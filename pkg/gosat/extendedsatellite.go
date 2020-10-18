package gosat

import (
	"math"
	"time"
)

type extendedSatellite struct {
	satellite
	Track [128][6]float64
}

// eci2ecef Conversion from Earth-centered inertial (ECI) to Earth-centered, Earth-Fixed (ECEF)
func eci2ecef(eci [6]float64, angle float64) (ecef [6]float64) {
	// angular speed of Earth's rotation rad/s
	angularSpeed := 7.2921151467e-5
	// alpha := angularSpeed*minutes()*60.0 + gsto()
	// alpha := 0.0
	sina := math.Sin(angle)
	cosa := math.Cos(angle)
	ecef[0] = cosa*eci[0] + sina*eci[1]
	ecef[1] = -sina*eci[0] + cosa*eci[1]
	ecef[2] = eci[2]
	ecef[3] = cosa*eci[3] + sina*eci[4] + angularSpeed*ecef[1]
	ecef[4] = -sina*eci[3] + cosa*eci[4] - angularSpeed*ecef[0]
	ecef[5] = eci[5]
	return
}

func (es *extendedSatellite) Update(t time.Time) error {
	timeStartSeconds := -60.0 * math.Pi / es.no
	timeStepSeconds := 120.0 * math.Pi / es.no / float64(len(es.Track))
	timeStart := time.Duration(float64(time.Second) * timeStartSeconds)
	timeStep := time.Duration(float64(time.Second) * timeStepSeconds)
	trackTime := t.Add(timeStart)
	for i := range es.Track {
		angle := float64(trackTime.Hour()*3600 + trackTime.Minute()*60 + trackTime.Second())
		angle += float64(trackTime.Nanosecond()) * 1.0e-9
		angle *= 7.2921151467e-5
		es.satellite.Update(trackTime)
		es.Track[i] = eci2ecef(es.Coords, angle)
		trackTime = trackTime.Add(timeStep)
	}
	// TODO make seconds of day func
	err := es.satellite.Update(t)
	angle := float64(t.Hour()*3600 + t.Minute()*60 + t.Second())
	angle += float64(t.Nanosecond()) * 1.0e-9
	angle *= 7.2921151467e-5
	es.Coords = eci2ecef(es.Coords, angle)
	return err
}
