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
			if tleItem.satn == 11801 {
				fmt.Print("Deep space, perigee = 82.48 (<98) for s4 > 20 mod")
			}
			err = sat.sgp4init("wgs84", tleItem, 'i')
			if err != nil {
				fmt.Println(*tleItem)
				t.Error(err)
			}
		}
	}
}
