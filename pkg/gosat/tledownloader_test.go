package gosat

import (
	"os"
	"testing"
)

func TestTleLoadList(t *testing.T) {
	err := loadList("testdata/tlelist.json")
	if err != nil {
		if os.IsNotExist(err) {
			t.Skip()
		} else {
			t.Error(err)
		}
	}
}
