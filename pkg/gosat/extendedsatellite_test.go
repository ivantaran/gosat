package gosat

import (
	"encoding/json"
	"testing"
)

func TestExtendedSatellite(t *testing.T) {
	var sat extendedSatellite
	tleList, err := LoadTle("testdata/SGP4-VER.TLE")
	if err != nil {
		t.Error(err)
	}
	err = sat.sgp4init("wgs84", &tleList[0], 'i')
	if err != nil {
		t.Fatal(err)
	}
	err = sat.Update(sat.Timestamp)
	if err != nil {
		t.Fatal(err)
	}
	bytes, err := json.Marshal(sat)
	if err != nil {
		t.Fatal(err)
	}
	t.Log(string(bytes))
}
