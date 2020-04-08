package gosat

import (
	"bufio"
	"fmt"
	"math"
	"os"
	"strconv"
	"strings"
)

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

func jday(year int, days float64) (jd float64, jdFrac float64) {
	mon, day, hr, minute, sec := days2mdhms(year, days)
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
func LoadTle(fileName string) ([]*Tle, error) {
	const xpdotp = 1440.0 / TwoPi

	file, err := os.Open(fileName)

	if err != nil {
		return nil, err
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)

	lineFirstOk := false
	lineSecondOk := false
	var list []*Tle
	t := new(Tle)
	for scanner.Scan() {
		line := scanner.Text()
		if len(line) <= 0 || line[0] == '#' {
			continue
		} else if len(line) >= 69 {
			switch line[0] {
			case '1':
				t.satn, err = strconv.Atoi(line[2:7])
				t.class = line[7]
				t.design = line[9:17]
				year, _ := strconv.Atoi(line[18:20])
				if year < 57 {
					year += 1900 // TODO aged units
				} else {
					year += 2000
				}
				days, _ := strconv.ParseFloat(line[20:32], 64)
				jd, jdFrac := jday(year, days)
				t.epoch = jd + jdFrac - 2433281.5
				t.xndot, err = strconv.ParseFloat(strings.TrimSpace(line[33:43]), 64)
				t.xnddot, err = strconv.ParseFloat(strings.TrimSpace(line[44:45]+"."+line[45:50]+"e"+line[50:52]), 64)
				t.xbstar, err = strconv.ParseFloat(strings.TrimSpace(line[53:54]+"."+line[54:59]+"e"+line[59:61]), 64)
				t.ephtype, err = strconv.Atoi(line[62:63])
				t.elsetn, err = strconv.Atoi(strings.TrimSpace(line[64:68]))
				t.cs1, err = strconv.Atoi(line[68:69])

				t.xno /= xpdotp
				t.xndot /= xpdotp * 1440.0
				t.xnddot /= xpdotp * 1440.0 * 1440.0

				lineFirstOk = true
				lineSecondOk = false
			case '2':
				satn, _ := strconv.Atoi(line[2:7])
				if satn != t.satn {
					fmt.Printf("different satn %d != %d\n", satn, t.satn)
				}
				t.xinclo, err = strconv.ParseFloat(strings.TrimSpace(line[8:16]), 64)
				t.xnodeo, err = strconv.ParseFloat(strings.TrimSpace(line[17:25]), 64)
				t.xecco, err = strconv.ParseFloat("."+strings.TrimSpace(line[26:33]), 64)
				t.xargpo, err = strconv.ParseFloat(strings.TrimSpace(line[34:42]), 64)
				t.xmo, err = strconv.ParseFloat(strings.TrimSpace(line[43:51]), 64)
				t.xno, err = strconv.ParseFloat(strings.TrimSpace(line[52:63]), 64)
				t.revn, err = strconv.Atoi(strings.TrimSpace(line[63:68]))
				t.cs2, err = strconv.Atoi(line[68:69])

				t.xinclo *= Deg2Rad
				t.xnodeo *= Deg2Rad
				t.xargpo *= Deg2Rad
				t.xmo *= Deg2Rad

				lineSecondOk = true
			default:
				break
			}
		} else {
			t.title = strings.TrimSpace(line)
		}
		if lineFirstOk && lineSecondOk {
			lineFirstOk = false
			lineSecondOk = false
			list = append(list, t)
			t = new(Tle)
		}
	}

	return list, nil
}
