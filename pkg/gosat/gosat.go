package gosat

import (
	"bufio"
	"encoding/json"
	"fmt"
	"io/ioutil"
	"log"
	"net"
	"os"
	"path/filepath"
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
	time         time.Time

	satMap map[string]*extendedSatellite
	tleMap map[string]Tle
	sun    sunObject
}

type inputStruct struct {
	IDList   []string
	AppendID []string
	RemoveID []string
	Time     time.Time
	Names    bool
}

//NewGosat TODO fill all fields here
func NewGosat() *Gosat {
	var gs Gosat
	gs.tleMap = make(map[string]Tle)
	gs.satMap = make(map[string]*extendedSatellite)
	gs.conns = make(map[*conn]struct{})
	return &gs
}

// loadIDList Load ID's List
func (gs *Gosat) loadIDList(bytes []byte, tle []Tle) error {
	var list inputStruct

	err := json.Unmarshal(bytes, &list)
	if err != nil {
		log.Fatal(err)
		return err
	}

	tleMap := make(map[string]Tle)
	for _, t := range tle {
		tleMap[t.Title] = t
	}

	for _, id := range list.IDList {
		if t, ok := tleMap[id]; ok {
			var sat extendedSatellite
			err := sat.sgp4init("wgs84", &t, 'i')
			if err != nil {
				log.Fatal(err)
				return err
			}
			gs.satMap[id] = &sat
		}
	}

	return nil
}

func (gs *Gosat) appendToSatMap(list []string) error {
	for _, id := range list {
		if t, ok := gs.tleMap[id]; ok {
			var sat extendedSatellite
			err := sat.sgp4init("wgs84", &t, 'i')
			if err != nil {
				log.Fatal(err)
				return err
			}
			gs.satMap[id] = &sat
		}
	}
	return nil
}

func (gs *Gosat) removeFromSatMap(list []string) error {
	for _, id := range list {
		if _, ok := gs.tleMap[id]; ok {
			delete(gs.satMap, id)
		}
	}
	return nil
}

func (gs *Gosat) update(time time.Time) (bytes []byte, err error) {
	for _, sat := range gs.satMap {
		err = sat.Update(time)
		if err != nil {
			fmt.Println(err)
		}
	}
	gs.sun.update(time)
	bytes, err = json.Marshal(struct {
		SatMap map[string]*extendedSatellite
		Sun    sunObject
	}{gs.satMap, gs.sun})
	return
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
	gs.conns[c] = struct{}{}
}

func (gs *Gosat) parseInput(bytes []byte) (packet []byte, err error) {
	iStruct := inputStruct{Names: false}
	err = json.Unmarshal(bytes, &iStruct)
	if len(iStruct.IDList) > 0 {
		gs.appendToSatMap(iStruct.IDList)
	}
	if len(iStruct.RemoveID) > 0 {
		gs.removeFromSatMap(iStruct.RemoveID)
	}
	if iStruct.Time.IsZero() {
		gs.time = time.Now().UTC()
	} else {
		gs.time = iStruct.Time
	}
	if iStruct.Names {
		packet, err = json.Marshal(struct{ TleMap map[string]Tle }{gs.tleMap})
		if err != nil {
			return
		}
	}
	return
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
	buffer := make([]byte, 0, 4096)
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
			buffer = append(buffer, scanr.Bytes()...)
			if json.Valid(buffer) {
				packet, err := gs.parseInput(buffer)
				if err != nil {
					log.Println(err)
					return err
				}
				if packet != nil {
					w.Write(packet)
					w.Write([]byte("\n"))
				}
				buffer = buffer[:0]

				bytes, err := gs.update(gs.time)
				if err != nil {
					return err
				}
				w.Write(bytes)
				w.Write([]byte("\n"))
				w.Flush()
			}
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

// LoadTle TODO comment me
func (gs *Gosat) LoadTle() (err error) {
	dir, err := os.UserCacheDir()
	if err != nil {
		return err
	}
	dir = filepath.Join(dir, "gosat")
	files, err := ioutil.ReadDir(dir)
	if err != nil {
		return err
	}
	for _, file := range files {
		log.Printf("Loading TLE: %s\n", file.Name())
		if file.IsDir() {
			break
		}
		filePath := filepath.Join(dir, file.Name())
		tleMap, err := LoadTleAsMap(filePath)
		if err != nil {
			return err
		}
		for id, tle := range tleMap {
			gs.tleMap[id] = tle
		}

	}
	return
}
