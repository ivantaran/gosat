package gosat

import (
	"bufio"
	"encoding/json"
	"io/ioutil"
	"net"
	"os"
	"testing"
	"time"
)

func TestLoadIDList(t *testing.T) {
	tleList, err := LoadTle("testdata/SGP4-VER.TLE")
	if err != nil {
		t.Error(err)
	}

	file, err := os.Open("testdata/idlist.json")
	if err != nil {
		t.Error(err)
	}
	defer file.Close()

	bytes, _ := ioutil.ReadAll(file)
	gs := NewGosat()
	err = gs.loadIDList(bytes, tleList)
	if err != nil {
		t.Error(err)
	}

	bytes, err = json.Marshal(gs.satMap)
	if err != nil {
		t.Error(err)
	}
	t.Log(string(bytes))
}

func TestServerProtectsAgaintSlowloris(t *testing.T) {
	gs := NewGosat()
	gs.IdleTimeout = 5 * time.Second
	gs.MaxReadBytes = 1000

	go gs.ListenAndServe()

	time.Sleep(1 * time.Second) // hack to wait for server to start
	conn, err := net.Dial("tcp", gs.Addr)
	if err != nil {
		t.Fatal(err)
	}
	// We slowly write to simulate a Slowloris. We should fail in one second
	// because we don't satisfy the application level requirement of sending a complete request (with newlines)
	// within 1 second
	w := bufio.NewWriter(conn)
	r := bufio.NewReader(conn)
	for {
		_, err = w.WriteString("{}\n")
		if err != nil {
			t.Fatal(err)
		}
		w.Flush()
		time.Sleep(200 * time.Millisecond)
		line, _, err := r.ReadLine()
		if err != nil {
			if err.Error() != "EOF" {
				t.Error(err)
			}
		} else if len(line) > 0 {
			t.Skip()
		}
	}
}

func TestMarshalSatellite(t *testing.T) {
	tleList, err := LoadTle("testdata/SGP4-VER.TLE")
	if err != nil {
		t.Error(err)
	}
	var sat satellite
	err = sat.sgp4init("wgs84", &tleList[0], 'i')
	bytes, err := json.Marshal(tleList)
	if err != nil {
		t.Fatal(err)
	}
	t.Log(string(bytes))
}
