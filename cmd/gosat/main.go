package main

import (
	"time"

	"../../pkg/gosat"
)

func main() {
	gs := gosat.Gosat{
		Addr:         ":8080",
		IdleTimeout:  5 * time.Second,
		MaxReadBytes: 1000,
	}
	gs.LoadTle()
	go gs.ListenAndServe()
	time.Sleep(10 * time.Second)
	select {}
}
