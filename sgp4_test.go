package gosat

import "testing"

func TestInit(t *testing.T) {
	t.Run("init test", func(t *testing.T) {
		var sat elsetrec
		var tle1 tle
		sat.sgp4init("wgs84", &tle1, 'i', 0, 0.0)
	})
}

func TestTleLoad(t *testing.T) {
	err := LoadTle("SGP4-VER.TLE")
	if err != nil {
		t.Error(err)
	}
}
