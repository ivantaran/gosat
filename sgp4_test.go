package gosat

import (
	"bufio"
	"fmt"
	"os"
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
			testMap[id] = make(coordLineMap)
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
	tleList, err := LoadTle("SGP4-VER.TLE")
	if err != nil {
		t.Error(err)
	} else {
		for _, tleItem := range tleList {
			if tleItem.satn == 33334 {
				fmt.Println("OLOLO")
			}
			err = sat.sgp4init("wgs84", tleItem, 'i')
			switch sat.satnum {
			case 33333:
				var r [6]float64
				err := sat.sgp4(3600.0, r[0:3], r[3:])
				if err == nil || sat.error != 4 {
					fmt.Println(*tleItem)
					t.Errorf("[FAIL] errnum:%d satnum:%d \"%s\"", sat.error, sat.satnum, err)
				} else {
					t.Logf("[ OK ] errnum:%d satnum:%d \"%s\"", sat.error, sat.satnum, err)
				}
			case 33334:
				if err == nil || sat.error != 3 {
					fmt.Println(*tleItem)
					t.Errorf("[FAIL] errnum:%d satnum:%d \"%s\"", sat.error, sat.satnum, err)
				} else {
					t.Logf("[ OK ] errnum:%d satnum:%d \"%s\"", sat.error, sat.satnum, err)
				}
			// case 33335: Wrong test case
			// 	var r [6]float64
			// 	err := sat.sgp4(3600.0, r[0:3], r[3:])
			// 	if err == nil || sat.error != 3 {
			// 		fmt.Println(*tleItem)
			// 		t.Errorf("[FAIL] errnum:%d satnum:%d \"%s\"", sat.error, sat.satnum, err)
			// 	} else {
			// 		t.Logf("[ OK ] errnum:%d satnum:%d \"%s\"", sat.error, sat.satnum, err)
			// 	}
			default:
				if err != nil {
					fmt.Println(*tleItem)
					t.Errorf("[FAIL] errnum:%d satnum:%d \"%s\"", sat.error, sat.satnum, err)
				}
			}
		}
	}
}

func TestReference(t *testing.T) {
	var r [6]float64

	testMap, err := loadRefTest("reftest.out")
	if err != nil {
		t.Error(err)
	} else {
		fmt.Println(len(testMap))
	}
	var sat elsetrec
	tleList, err := LoadTle("SGP4-VER.TLE")
	if err != nil {
		t.Error(err)
	} else {
		for _, tleItem := range tleList {
			err = sat.sgp4init("wgs84", tleItem, 'i')
			if coordLineMap, ok := testMap[sat.satnum]; ok {
				for time, coord := range coordLineMap {
					sat.sgp4(time, r[0:3], r[3:])
					fmt.Printf("%d %f %f %f\n", sat.satnum, time, coord[0], r[0])
				}
			}
		}
	}
}
