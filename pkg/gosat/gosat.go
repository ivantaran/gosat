package gosat

import (
	"bufio"
	"encoding/json"
	"log"
	"net"
	"strings"
	"sync"
	"time"
)

// Gosat TODO comment me
type Gosat struct {
	Addr         string
	conns        map[*conn]struct{}
	IdleTimeout  time.Duration
	inShutdown   bool
	listener     net.Listener
	mu           sync.Mutex
	MaxReadBytes int64

	satMap map[int]elsetrec
}

type idList struct {
	IDList []int `json:"idList"`
}

// loadIDList Load ID's List
func (gs *Gosat) loadIDList(bytes []byte, tle []Tle) error {
	var list idList

	err := json.Unmarshal(bytes, &list)
	if err != nil {
		log.Fatal(err)
		return err
	}

	tleMap := make(map[int]Tle)
	for _, t := range tle {
		id := t.Satnum
		tleMap[id] = t
	}

	gs.satMap = make(map[int]elsetrec)
	for _, id := range list.IDList {
		if t, ok := tleMap[id]; ok {
			var sat elsetrec
			err := sat.sgp4init("wgs84", &t, 'i')
			if err != nil {
				log.Fatal(err)
				return err
			}
			gs.satMap[id] = sat
		}
	}

	return nil
}

func (gs *Gosat) update(t float64) {
	for _, sat := range gs.satMap {
		sat.sgp4(t)
	}
}

// ListenAndServe TODO comment me
// link to reference design: https://github.com/sahilm/shouter/
func (gs *Gosat) ListenAndServe() error {
	addr := gs.Addr
	if addr == "" {
		addr = ":8080"
	}
	log.Printf("starting server on %v\n", addr)
	listener, err := net.Listen("tcp", addr)
	if err != nil {
		return err
	}
	defer listener.Close()
	gs.listener = listener
	for {
		// should be guarded by mu
		if gs.inShutdown {
			break
		}
		newConn, err := listener.Accept()
		if err != nil {
			log.Printf("error accepting connection %v", err)
			continue
		}
		log.Printf("accepted connection from %v", newConn.RemoteAddr())
		conn := &conn{
			Conn:          newConn,
			IdleTimeout:   gs.IdleTimeout,
			MaxReadBuffer: gs.MaxReadBytes,
		}
		gs.trackConn(conn)
		conn.SetDeadline(time.Now().Add(conn.IdleTimeout))
		go gs.handle(conn)
	}
	return nil
}

func (gs *Gosat) trackConn(c *conn) {
	defer gs.mu.Unlock()
	gs.mu.Lock()
	if gs.conns == nil {
		gs.conns = make(map[*conn]struct{})
	}
	gs.conns[c] = struct{}{}
}

func (gs *Gosat) handle(c *conn) error {
	defer func() {
		log.Printf("closing connection from %v", c.RemoteAddr())
		c.Close()
		gs.deleteConn(c)
	}()
	r := bufio.NewReader(c)
	w := bufio.NewWriter(c)
	scanr := bufio.NewScanner(r)

	sc := make(chan bool)
	deadline := time.After(c.IdleTimeout)
	for {
		go func(s chan bool) {
			s <- scanr.Scan()
		}(sc)
		select {
		case <-deadline:
			return nil
		case scanned := <-sc:
			if !scanned {
				if err := scanr.Err(); err != nil {
					return err
				}
				return nil
			}
			w.WriteString(strings.ToUpper(scanr.Text()) + "\n")
			w.Flush()
			deadline = time.After(c.IdleTimeout)
		}
	}
}

func (gs *Gosat) deleteConn(c *conn) {
	defer gs.mu.Unlock()
	gs.mu.Lock()
	delete(gs.conns, c)
}

// Shutdown TODO comment me
func (gs *Gosat) Shutdown() {
	// should be guarded by mu
	gs.inShutdown = true
	log.Println("shutting down...")
	gs.listener.Close()
	ticker := time.NewTicker(500 * time.Millisecond)
	defer ticker.Stop()
	for {
		select {
		case <-ticker.C:
			log.Printf("waiting on %v connections", len(gs.conns))
		}
		if len(gs.conns) == 0 {
			return
		}
	}
}
