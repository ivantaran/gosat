package gosat

import (
	"bufio"
	"fmt"
	"math"
	"os"
	"sort"
	"strconv"
	"strings"
	"testing"
)

type coordLine []float64
type coordLineMap map[float64]coordLine
type testMap map[int]coordLineMap

func loadRefTest(fileName string) (testMap, error) {
	file, err := os.Open(fileName)
	if err != nil {
		return nil, err
	}
	defer file.Close()
	scanner := bufio.NewScanner(file)
	testMap := make(testMap)

	id := 0
	for scanner.Scan() {
		line := scanner.Text()
		fields := strings.Fields(line)
		if len(fields) == 2 && fields[1] == "xx" {
			id, err = strconv.Atoi(fields[0])
			if err != nil {
				return nil, err
			}
			if testMap[id] == nil {
				testMap[id] = make(coordLineMap)
			}
		} else if len(fields) > 20 {
			fmt.Println("last long line test")
		} else if len(fields) >= 7 {
			time, err := strconv.ParseFloat(fields[0], 64)
			if err != nil {
				return nil, err
			}
			data := make([]float64, 6)
			for i, field := range fields[1:] {
				data[i], err = strconv.ParseFloat(field, 64)
				if err != nil {
					return nil, err
				}
				if i == 5 {
					break
				}
			}
			testMap[id][time] = data
		}
	}

	return testMap, nil
}

func TestInit(t *testing.T) {
	t.Run("init test", func(t *testing.T) {
		var sat elsetrec
		var tle1 Tle
		sat.sgp4init("wgs84", &tle1, 'i')
	})
}

func TestTleLoad(t *testing.T) {
	var sat elsetrec
	tleList, err := LoadTle("./testdata/SGP4-VER.TLE")
	if err != nil {
		t.Error(err)
	}
	for _, tleItem := range tleList {
		err = sat.sgp4init("wgs84", &tleItem, 'i')
		switch sat.Satnum {
		case 33333:
			err := sat.sgp4(3600.0)
			if err == nil || sat.error != 4 {
				fmt.Println(tleItem)
				t.Errorf("[FAIL] errnum:%d satnum:%d \"%s\"", sat.error, sat.Satnum, err)
			} else {
				t.Logf("[ OK ] errnum:%d satnum:%d \"%s\"", sat.error, sat.Satnum, err)
			}
		case 33334:
			if err == nil || sat.error != 3 {
				fmt.Println(tleItem)
				t.Errorf("[FAIL] errnum:%d satnum:%d \"%s\"", sat.error, sat.Satnum, err)
			} else {
				t.Logf("[ OK ] errnum:%d satnum:%d \"%s\"", sat.error, sat.Satnum, err)
			}
		// case 33335: Wrong test case
		// 	var r [6]float64
		// 	err := sat.sgp4(3600.0)
		// 	if err == nil || sat.error != 3 {
		// 		fmt.Println(*tleItem)
		// 		t.Errorf("[FAIL] errnum:%d satnum:%d \"%s\"", sat.error, sat.Satnum, err)
		// 	} else {
		// 		t.Logf("[ OK ] errnum:%d satnum:%d \"%s\"", sat.error, sat.Satnum, err)
		// 	}
		default:
			if err != nil {
				fmt.Println(tleItem)
				t.Errorf("[FAIL] errnum:%d satnum:%d \"%s\"", sat.error, sat.Satnum, err)
			}
		}
	}
}

func TestReference(t *testing.T) {
	gravConst := "wgs72"
	treshold := 1.0e-8
	testFileName := "./testdata/reftest_java.out"
	// testFileName := "reftest_cpp.out"
	testMap, err := loadRefTest(testFileName)
	if err != nil {
		t.Error(err)
	} else {
		if len(testMap) != 32 {
			t.Errorf("len(testMap): %d != 32\n", len(testMap))
		}
	}
	var sat elsetrec
	tleList, err := LoadTle("./testdata/SGP4-VER.TLE")
	if err != nil {
		t.Error(err)
	}
	for _, tleItem := range tleList {
		// err = sat.sgp4init("wgs84", tleItem, 'i')
		err = sat.sgp4init(gravConst, &tleItem, 'i')
		if coordLineMap, ok := testMap[sat.Satnum]; ok {
			times := make([]float64, 0.0, len(coordLineMap))
			for time := range coordLineMap {
				times = append(times, time)
			}
			sort.Float64s(times)
			for _, time := range times {
				coord := coordLineMap[time]
				sat.sgp4(time)
				ok := true
				for i := 0; i < len(sat.Coords); i++ {
					delta := coord[i] - sat.Coords[i]
					if math.Abs(delta) > treshold {
						ok = false
						break
					}
				}
				if !ok {
					fmt.Printf("satnum:%d time:%0.1f\n", sat.Satnum, time)
					for i := 0; i < len(sat.Coords); i++ {
						fmt.Printf("\t%+0.4e %+0.4e %+0.4e\n", coord[i], sat.Coords[i], coord[i]-sat.Coords[i])
					}
					t.Error()
				}
			}
		}
	}
}
