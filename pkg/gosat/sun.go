package gosat

import (
	"math"
	"time"
)

const (
	degToRad = math.Pi / 180.0
	radToDeg = 180.0 / math.Pi
)

func timeJulianCent(jd float64) float64 {
	return (jd - 2451545.0) / 36525.0
}

func jdFromJulianCent(t float64) float64 {
	return t*36525.0 + 2451545.0
}

func isLeapYear(yr int) bool {
	return (yr%4 == 0 && yr%100 != 0) || (yr%400 == 0)
}

func doyFromJd(jd float64) float64 {

	z := math.Floor(jd + 0.5)
	f := (jd + 0.5) - z

	a := z
	if z >= 2299161.0 {
		alpha := math.Floor((z - 1867216.25) / 36524.25)
		a += 1.0 + alpha - math.Floor(alpha*0.25)
	}

	b := a + 1524.0
	c := math.Floor((b - 122.1) / 365.25)
	d := math.Floor(365.25 * c)
	e := math.Floor((b - d) / 30.6001)

	day := b - d - math.Floor(30.6001*e) + f
	month := e
	if e < 14.0 {
		month -= 1.0
	} else {
		month -= 13.0
	}
	year := c
	if month > 2.0 {
		year -= 4716.0
	} else {
		year -= 4715.0
	}

	var k float64
	if isLeapYear(int(year)) {
		k = 1.0
	} else {
		k = 2.0
	}

	return math.Floor((275.0*month)/9.0) - k*math.Floor((month+9.0)/12.0) + day - 30.0
}

func geomMeanLongSun(t float64) float64 {
	lon := 280.46646 + t*(36000.76983+t*0.0003032)
	lon = math.Mod(lon, 360.0)
	if lon < 0.0 {
		lon += 360.0
	}
	lon = lon * degToRad
	return lon // in radians
}

func geom_mean_anomaly_sun(t float64) float64 {
	m := degToRad * (357.52911 + t*(35999.05029-0.0001537*t))
	return m // in radians
}

func eccentricity_earth_orbit(t float64) float64 {
	e := 0.016708634 - t*(0.000042037+0.0000001267*t)
	return e // unitless
}

func sun_rad_vector(t, m, seoc, e float64) float64 {
	v := m + seoc
	r := (1.000001018 * (1.0 - e*e)) / (1.0 + e*math.Cos(v))
	return r // in AUs
}

func sun_apparent_long(t, true_long float64) float64 {
	omega := degToRad * (125.04 - 1934.136*t)
	lambda := true_long - degToRad*(0.00569+0.00478*math.Sin(omega))
	return lambda // in radians
}

func mean_obliquity_of_ecliptic(t float64) float64 {
	seconds := 21.448 - t*(46.8150+t*(0.00059-t*0.001813))
	e0 := degToRad * (23.0 + (26.0+seconds/60.0)/60.0)
	return e0 // in radians
}

func obliquity_correction(t float64) float64 {
	e0 := mean_obliquity_of_ecliptic(t)
	omega := degToRad * (125.04 - 1934.136*t)
	e := e0 + degToRad*(0.00256*math.Cos(omega))
	return e // in radians
}

func sun_declination(t, true_long float64) float64 {
	e := obliquity_correction(t)
	lambda := sun_apparent_long(t, true_long)

	sint := math.Sin(e) * math.Sin(lambda)
	theta := math.Asin(sint)
	return theta // in radians
}

func get_jd( /*struct tm *gmt*/ ) float64 {
	// int month, year;
	// float64 a, b, jd;

	// month = gmt->tm_mon + 1;
	// year = gmt->tm_year + 1900;

	// if (month <= 2) {
	//     year -= 1;
	//     month += 12;
	// }

	// a = floor(year * 0.01);
	// b = 2.0 - a + floor(a * 0.25);
	// jd = floor(365.25 * (year + 4716.0)) + floor(30.6001 * (month + 1.0))
	//         + gmt->tm_mday + b - 1524.5;
	// return jd;
	return 0.0
}

func atmo_refraction_correction(elevation float64) float64 {
	var value, te float64
	/* Atmospheric refraction correction */
	if elevation > degToRad*(85.0) {
		value = 0.0
	} else {
		te = math.Tan(elevation)
		if elevation > degToRad*(5.0) {
			value = 58.1/te - 0.07/(te*te*te) +
				0.000086/(te*te*te*te*te)
		} else {
			if elevation > degToRad*(-0.575) {
				value = 1735.0 + elevation*(-518.2+elevation*
					(103.4+elevation*(-12.79+elevation*0.711)))
			} else {
				value = -20.774 / te
			}
		}
		value = degToRad * (value / 3600.0)
	}
	return value
}

// ae0 useArc - atmospheric refraction correction flag
func ae0(slat, slon, lat, lon float64, useArc rune) (azimuth, elevation float64) {

	omega := lon - slon
	sinlat := math.Sin(lat)
	coslat := math.Cos(lat)
	sinslat := math.Sin(slat)
	cosslat := math.Cos(slat)

	csz := sinlat*sinslat + coslat*cosslat*math.Cos(omega)

	if csz > 1.0 {
		csz = 1.0
	} else {
		if csz < -1.0 {
			csz = -1.0
		}
	}

	elevation = math.Asin(csz)
	azDenom := coslat * math.Cos(elevation)
	azimuth = (sinlat*csz - sinslat) / azDenom

	if math.IsInf(azimuth, 0) {
		if lat > 0.0 {
			azimuth = math.Pi
		} else {
			azimuth = 0.0
		}
	} else {
		if math.Abs(azimuth) > 1.0 {
			if azimuth < 0.0 {
				azimuth = -1.0
			} else {
				azimuth = 1.0
			}
		}

		azimuth = math.Pi - math.Acos(azimuth)
		if omega > 0.0 {
			azimuth = -azimuth
		}
	}

	if azimuth < 0.0 {
		azimuth += 2.0 * math.Pi
	}

	if useArc == 1 || useArc == 'y' {
		elevation += atmo_refraction_correction(elevation)
	}

	//    if (elevation < degToRad * (-18.0)) {
	//        puts("A Night at the Roxbury");
	//    }

	return
}

func ll0(t, localtime float64) (latitude, longitude float64) {

	m := geom_mean_anomaly_sun(t)
	e := eccentricity_earth_orbit(t)

	// equation_of_time
	epsilon := obliquity_correction(t)
	l0 := geomMeanLongSun(t)

	y := math.Tan(epsilon * 0.5)
	y *= y

	sin2l0 := math.Sin(2.0 * l0)
	sin4l0 := math.Sin(4.0 * l0)
	cos2l0 := math.Cos(2.0 * l0)
	sinm := math.Sin(m)
	sin2m := math.Sin(2.0 * m)
	sin3m := math.Sin(3.0 * m)

	etime := y*sin2l0 - 2.0*e*sinm + 4.0*e*y*sinm*cos2l0 - 0.5*y*y*sin4l0 - 1.25*e*e*sin2m
	eqTime := radToDeg * (etime) * 4.0 // in minutes of time
	// equation_of_time

	// sun_eq_of_center
	seoc := sinm*(1.914602-t*(0.004817+0.000014*t)) + sin2m*(0.019993-0.000101*t) + sin3m*0.000289
	seoc = degToRad * (seoc)
	// sun_eq_of_center

	trueLong := seoc + l0

	latitude = sun_declination(t, trueLong)
	longitude = 180.0 - (eqTime+localtime)*0.25
	longitude = math.Mod(longitude+180.0, 360.0)
	if longitude < 0.0 {
		longitude += 360.0
	}
	longitude -= 180.0
	longitude = degToRad * (longitude)

	return
}

func ll(t time.Time) (latitude, longitude float64) {
	jday := julian(t)
	julianCenturies := timeJulianCent(jday)
	minutesOfTheDay := secondsOfTheDay(t) / 60.0
	// TODO julianCenturies must be at 12:00
	latitude, longitude = ll0(julianCenturies, minutesOfTheDay)
	return
}

func ae(t time.Time, latitude, longitude float64, useArc rune) (azimuth, elevation float64) {
	sunLatitude, sunLongitude := ll(t)
	azimuth, elevation = ae0(sunLatitude, sunLongitude, latitude, longitude, useArc)
	return
}
