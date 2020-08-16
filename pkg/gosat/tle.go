package gosat

import (
	"bufio"
	"fmt"
	"math"
	"os"
	"strconv"
	"strings"
	"time"
)

//Tle TLE record container
type Tle struct {
	/* Auxiliary fields */
	class     byte   // Classification (U=Unclassified, C=Classified, S=Secret)
	design    string // International Designator YYNNNAAA (YY - launch year, NNN - launch number of the year, AAA - piece of the launch)
	elsetn    int    // Element set number. Incremented when a new TLE is generated for this object
	ephtype   int    // Ephemeris type (internal use only - always zero in distributed TLE data)
	revnum    int    // Revolution number at epoch
	Satnum    int    // Satellite catalog number
	Title     string // Title line, length 24
	Timestamp time.Time
	/* Orbit parameters */
	argpo float64 // Argument of Perigee
	bstar float64 // Drag Term aka Radiation Pressure Coefficient or BSTAR
	ecco  float64 // Eccentricity
	epoch float64 // Epoch ???
	inclo float64 // Inclination
	mo    float64 // Mean Anomaly
	nddot float64 // Second Derivative of Mean Motion
	ndot  float64 // First Derivative of Mean Motion aka the Ballistic Coefficient
	no    float64 // Mean Motion
	nodeo float64 // Right Ascension of the Ascending Node
}

func days2mdhms(year int, days float64) (
	mon int, day int, hr int, minute int, sec float64) {

	var lmonth = []int{0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31}
	dayofyr := int(days)

	/* ----------------- find month and day of month ---------------- */
	if (year % 4) == 0 {
		lmonth[2] = 29
	}

	inttemp := int(0)
	i := 0
	for (dayofyr > inttemp+lmonth[i]) && (i < 12) {
		inttemp = inttemp + lmonth[i]
		i++
	}

	mon = i
	day = dayofyr - inttemp

	/* ----------------- find hours minutes and seconds ------------- */
	temp := (days - float64(dayofyr)) * 24.0
	hr = int(temp)
	temp = (temp - float64(hr)) * 60.0
	minute = int(temp)
	sec = (temp - float64(minute)) * 60.0

	return
}

func jday(year int, days float64) (jd float64, jdFrac float64, timestamp time.Time) {
	mon, day, hr, minute, sec := days2mdhms(year, days)
	isec, fsec := math.Modf(sec)
	timestamp = time.Date(year, time.Month(mon), day, hr, minute, int(isec), int(fsec*1.0e9), time.UTC)
	jd = 367.0*float64(year) -
		math.Floor((7.0*(float64(year)+math.Floor(float64(mon+9)/12.0)))*0.25) +
		math.Floor(275.0*float64(mon)/9.0) +
		float64(day) + 1721013.5 // use - 678987.0 to go to mjd directly
	jdFrac = (sec + float64(minute*60+hr*3600)) / 86400.0

	// check that the day and fractional day are correct
	if math.Abs(jdFrac) > 1.0 {
		dtt := math.Floor(jdFrac)
		jd = jd + dtt
		jdFrac = jdFrac - dtt
	}
	// - 0.5*sgn(100.0*year + mon - 190002.5) + 0.5;
	return
}

// LoadTle load TLE data
func LoadTle(fileName string) ([]Tle, error) {
	const xpdotp = 1440.0 / twoPi

	file, err := os.Open(fileName)

	if err != nil {
		return nil, err
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)

	lineFirstOk := false
	lineSecondOk := false
	var list []Tle
	t := Tle{}
	for scanner.Scan() {
		line := scanner.Text()
		if len(line) <= 0 || line[0] == '#' {
			continue
		} else if len(line) >= 69 {
			switch line[0] {
			case '1':
				t.Satnum, err = strconv.Atoi(line[2:7])
				t.class = line[7]
				t.design = line[9:17]
				year, _ := strconv.Atoi(line[18:20])
				if year < 57 {
					year += 2000 // TODO aged units
				} else {
					year += 1900
				}
				days, _ := strconv.ParseFloat(line[20:32], 64)
				jd, jdFrac, timestamp := jday(year, days)
				t.Timestamp = timestamp
				t.epoch = jd + jdFrac - 2433281.5
				t.ndot, err = strconv.ParseFloat(strings.TrimSpace(line[33:43]), 64)
				t.nddot, err = strconv.ParseFloat(strings.TrimSpace(line[44:45]+"."+line[45:50]+"e"+line[50:52]), 64)
				t.bstar, err = strconv.ParseFloat(strings.TrimSpace(line[53:54]+"."+line[54:59]+"e"+line[59:61]), 64)
				t.ephtype, err = strconv.Atoi(line[62:63])
				t.elsetn, err = strconv.Atoi(strings.TrimSpace(line[64:68]))
				// cs1, err = strconv.Atoi(line[68:69])

				t.ndot /= xpdotp * 1440.0
				t.nddot /= xpdotp * 1440.0 * 1440.0

				lineFirstOk = true
				lineSecondOk = false
			case '2':
				satn, _ := strconv.Atoi(line[2:7])
				if satn != t.Satnum {
					fmt.Printf("different satn %d != %d\n", satn, t.Satnum)
				}
				t.inclo, err = strconv.ParseFloat(strings.TrimSpace(line[8:16]), 64)
				t.nodeo, err = strconv.ParseFloat(strings.TrimSpace(line[17:25]), 64)
				t.ecco, err = strconv.ParseFloat("."+strings.TrimSpace(line[26:33]), 64)
				t.argpo, err = strconv.ParseFloat(strings.TrimSpace(line[34:42]), 64)
				t.mo, err = strconv.ParseFloat(strings.TrimSpace(line[43:51]), 64)
				t.no, err = strconv.ParseFloat(strings.TrimSpace(line[52:63]), 64)
				t.revnum, err = strconv.Atoi(strings.TrimSpace(line[63:68]))
				// cs2, err = strconv.Atoi(line[68:69])

				t.inclo *= deg2rad
				t.nodeo *= deg2rad
				t.argpo *= deg2rad
				t.mo *= deg2rad
				t.no /= xpdotp

				lineSecondOk = true
			default:
				break
			}
		} else {
			t.Title = strings.TrimSpace(line)
		}
		if lineFirstOk && lineSecondOk {
			lineFirstOk = false
			lineSecondOk = false
			list = append(list, t)
			// t := Tle{}
		}
	}

	return list, nil
}

// LoadTleAsMap load TLE data
func LoadTleAsMap(fileName string) (map[int]Tle, error) {
	const xpdotp = 1440.0 / twoPi

	file, err := os.Open(fileName)

	if err != nil {
		return nil, err
	}
	defer file.Close()

	tleMap := make(map[int]Tle)
	scanner := bufio.NewScanner(file)

	lineFirstOk := false
	lineSecondOk := false

	t := Tle{}
	for scanner.Scan() {
		line := scanner.Text()
		if len(line) <= 0 || line[0] == '#' {
			continue
		} else if len(line) >= 69 {
			switch line[0] {
			case '1':
				t.Satnum, err = strconv.Atoi(line[2:7])
				t.class = line[7]
				t.design = line[9:17]
				year, _ := strconv.Atoi(line[18:20])
				if year < 57 {
					year += 2000 // TODO aged units
				} else {
					year += 1900
				}
				days, _ := strconv.ParseFloat(line[20:32], 64)
				jd, jdFrac, timestamp := jday(year, days)
				t.Timestamp = timestamp
				t.epoch = jd + jdFrac - 2433281.5
				t.ndot, err = strconv.ParseFloat(strings.TrimSpace(line[33:43]), 64)
				t.nddot, err = strconv.ParseFloat(strings.TrimSpace(line[44:45]+"."+line[45:50]+"e"+line[50:52]), 64)
				t.bstar, err = strconv.ParseFloat(strings.TrimSpace(line[53:54]+"."+line[54:59]+"e"+line[59:61]), 64)
				t.ephtype, err = strconv.Atoi(line[62:63])
				t.elsetn, err = strconv.Atoi(strings.TrimSpace(line[64:68]))
				// cs1, err = strconv.Atoi(line[68:69])

				t.ndot /= xpdotp * 1440.0
				t.nddot /= xpdotp * 1440.0 * 1440.0

				lineFirstOk = true
				lineSecondOk = false
			case '2':
				satn, _ := strconv.Atoi(line[2:7])
				if satn != t.Satnum {
					fmt.Printf("different satn %d != %d\n", satn, t.Satnum)
				}
				t.inclo, err = strconv.ParseFloat(strings.TrimSpace(line[8:16]), 64)
				t.nodeo, err = strconv.ParseFloat(strings.TrimSpace(line[17:25]), 64)
				t.ecco, err = strconv.ParseFloat("."+strings.TrimSpace(line[26:33]), 64)
				t.argpo, err = strconv.ParseFloat(strings.TrimSpace(line[34:42]), 64)
				t.mo, err = strconv.ParseFloat(strings.TrimSpace(line[43:51]), 64)
				t.no, err = strconv.ParseFloat(strings.TrimSpace(line[52:63]), 64)
				t.revnum, err = strconv.Atoi(strings.TrimSpace(line[63:68]))
				// cs2, err = strconv.Atoi(line[68:69])

				t.inclo *= deg2rad
				t.nodeo *= deg2rad
				t.argpo *= deg2rad
				t.mo *= deg2rad
				t.no /= xpdotp

				lineSecondOk = true
			default:
				break
			}
		} else {
			t.Title = strings.TrimSpace(line)
		}
		if lineFirstOk && lineSecondOk {
			lineFirstOk = false
			lineSecondOk = false
			tleMap[t.Satnum] = t
			// t := Tle{}
		}
	}

	return tleMap, nil
}

func (dst *Tle) copy(src *Tle) {
	dst.class = src.class
	dst.design = src.design
	dst.elsetn = src.elsetn
	dst.ephtype = src.ephtype
	dst.revnum = src.revnum
	dst.Satnum = src.Satnum
	dst.Title = src.Title
	dst.Timestamp = src.Timestamp

	dst.argpo = src.argpo
	dst.bstar = src.bstar
	dst.ecco = src.ecco
	dst.epoch = src.epoch
	dst.inclo = src.inclo
	dst.mo = src.mo
	dst.nddot = src.nddot
	dst.ndot = src.ndot
	dst.no = src.no
	dst.nodeo = src.nodeo
}
