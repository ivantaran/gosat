package gosat

import (
	"io/ioutil"
	"os"
	"testing"
)

func TestLoadIDList(t *testing.T) {
	tleList, err := LoadTle("testdata/SGP4-VER.TLE")
	if err != nil {
		t.Error(err)
	}

	file, err := os.Open("testdata/idlist.json")
	if err != nil {
		t.Error(err)
	}
	defer file.Close()

	bytes, _ := ioutil.ReadAll(file)

	err = LoadIDList(bytes, tleList)
	if err != nil {
		t.Error(err)
	}
}
