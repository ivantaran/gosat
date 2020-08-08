package gosat

import (
	"encoding/json"
	"log"
)

// Gosat TODO comment me
type Gosat struct {
	satMap map[int]elsetrec
}

type idList struct {
	IDList []int `json:"idList"`
}

// loadIDList Load ID's List
func (gs *Gosat) loadIDList(bytes []byte, tle []Tle) error {
	var list idList

	err := json.Unmarshal(bytes, &list)
	if err != nil {
		log.Fatal(err)
		return err
	}

	tleMap := make(map[int]Tle)
	for _, t := range tle {
		id := t.satnum
		tleMap[id] = t
	}

	gs.satMap = make(map[int]elsetrec)
	for _, id := range list.IDList {
		if t, ok := tleMap[id]; ok {
			var sat elsetrec
			err := sat.sgp4init("wgs84", &t, 'i')
			if err != nil {
				log.Fatal(err)
				return err
			}
			gs.satMap[id] = sat
		}
	}

	return nil
}

func (gs *Gosat) update(t float64) {
	r := []float64{0.0, 0.0, 0.0}
	v := []float64{0.0, 0.0, 0.0}
	for _, sat := range gs.satMap {
		sat.sgp4(t, r, v)
	}
}
