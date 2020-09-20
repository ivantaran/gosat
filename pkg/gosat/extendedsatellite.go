package gosat

import (
	"math"
	"time"
)

type extendedSatellite struct {
	satellite
	Track [128][6]float64
}

func (es *extendedSatellite) Update(t time.Time) error {
	timeStartSeconds := -60.0 * math.Pi / es.no
	timeStepSeconds := 120.0 * math.Pi / es.no / float64(len(es.Track))
	timeStart := time.Duration(float64(time.Second) * timeStartSeconds)
	timeStep := time.Duration(float64(time.Second) * timeStepSeconds)
	trackTime := t.Add(timeStart)
	for i := range es.Track {
		es.satellite.Update(trackTime)
		es.Track[i] = es.Coords
		trackTime = trackTime.Add(timeStep)
	}
	err := es.satellite.Update(t)
	return err
}
