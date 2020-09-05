package main

import (
	"time"

	"../../pkg/gosat"
)

func main() {
	gs := gosat.NewGosat()
	gs.Addr = ":8080"
	gs.IdleTimeout = 30 * time.Second
	gs.MaxReadBytes = 1000
	gs.LoadTle()
	go gs.ListenAndServe()
	time.Sleep(10 * time.Second)
	select {}
}
