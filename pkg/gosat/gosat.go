package gosat

import (
	"encoding/json"
	"log"
)

type idList struct {
	IDList []int `json:"idList"`
}

// LoadIDList Load ID's List
func LoadIDList(bytes []byte, tle []*Tle) error {
	var list idList

	err := json.Unmarshal(bytes, &list)
	if err != nil {
		log.Fatal(err)
		return err
	}

	tleMap := make(map[int]*Tle)
	for _, t := range tle {
		id := t.satnum
		tleMap[id] = t
	}

	satMap := make(map[int]elsetrec)
	for _, id := range list.IDList {
		if t, ok := tleMap[id]; ok {
			var sat elsetrec
			err := sat.sgp4init("wgs84", t, 'i')
			if err != nil {
				log.Fatal(err)
				return err
			}
			satMap[id] = sat
		}
	}

	return nil
}
