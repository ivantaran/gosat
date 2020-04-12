package gosat

import (
	"fmt"
	"testing"
)

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
