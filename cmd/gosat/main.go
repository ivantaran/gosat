package main

import (
	"time"

	"../../pkg/gosat"
)

func main() {
	gs := gosat.NewGosat()
	gs.LoadTle()
	go gs.ListenAndServe()
	time.Sleep(10 * time.Second)
	select {}
}
