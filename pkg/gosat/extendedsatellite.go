package gosat

import (
	"math"
	"time"
)

const (
	// angular speed of Earth's rotation rad/s
	angularSpeed = 7.2921151467e-5

	// reference time for Julian Date 01.01.2020 12:00:00 UTC
	// http://www.onlineconversion.com/julian_date.htm
	julianDate2020 = 2458850.0
)

var (
	// reference time for Julian Date 01.01.2020 12:00:00 UTC
	unixJulianDate2020 = time.Unix(1577880000, 0)
	oneDay             = float64(86400.0 * time.Second)
)

type extendedSatellite struct {
	satellite
	Track [128][6]float64
}

// eci2ecef Conversion from Earth-centered inertial (ECI) to Earth-centered, Earth-Fixed (ECEF)
func eci2ecef(eci [6]float64, angle float64) (ecef [6]float64) {
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

func julian(t time.Time) float64 {
	return julianDate2020 + float64(t.Sub(unixJulianDate2020))/oneDay
}

func secondsOfTheDay(t time.Time) float64 {
	seconds := float64(t.Hour()*3600 + t.Minute()*60 + t.Second())
	seconds += float64(t.Nanosecond()) * 1.0e-9
	return seconds
}

func (es *extendedSatellite) Update(t time.Time) error {
	timeStartSeconds := -60.0 * math.Pi / es.no
	timeStepSeconds := 120.0 * math.Pi / es.no / float64(len(es.Track))
	timeStart := time.Duration(float64(time.Second) * timeStartSeconds)
	timeStep := time.Duration(float64(time.Second) * timeStepSeconds)
	trackTime := t.Add(timeStart)
	for i := range es.Track {
		jd := julian(trackTime)
		gst := gstime(jd)
		es.satellite.Update(trackTime)
		es.Track[i] = eci2ecef(es.Coords, gst)
		trackTime = trackTime.Add(timeStep)
	}
	// if 'julian' returned julian date without fractional part of the day
	// then rotation angle would be:
	// angle := secondsOfTheDay(trackTime)*angularSpeed + gst
	jd := julian(t)
	gst := gstime(jd)
	err := es.satellite.Update(t)
	es.Coords = eci2ecef(es.Coords, gst)
	return err
}
