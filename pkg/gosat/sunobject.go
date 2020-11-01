package gosat

import "time"

type sunObject struct {
	Coords [6]float64
}

func (so *sunObject) update(time time.Time) error {
	so.Coords = SunEcef(time)
	return nil
}
