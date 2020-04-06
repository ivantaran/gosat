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
func LoadTle(fileName string) error {

	file, err := os.Open(fileName)

	if err != nil {
		return err
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	var t tle
	for i := 0; i < 10 && scanner.Scan(); i++ {
		line := scanner.Text()
		if len(line) <= 0 {
			continue
		}
		switch line[0] {
		case '#':
			fmt.Println("comment")
		case '1':
			fmt.Println("line 1")
			t.satn, _ = strconv.Atoi(line[2:7])
			t.class = line[7]
			t.design = line[9:17]
			year, _ := strconv.Atoi(line[18:20])
			days, _ := strconv.ParseFloat(line[20:32], 64)
			jd, jdFrac := jday(year, days)
			t.epoch = jd + jdFrac
			t.xndot, _ = strconv.ParseFloat(line[33:43], 64)
			t.xnddot, _ = strconv.ParseFloat(line[44:45]+"."+line[45:50]+"e"+line[50:52], 64)
			fmt.Println(t, string(t.class))
		case '2':
			fmt.Println("line 2")
		default:
			fmt.Println(strings.TrimSpace(line))
		}
	}

	return nil
}
