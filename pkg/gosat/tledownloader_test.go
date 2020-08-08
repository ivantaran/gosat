package gosat

import (
	"fmt"
	"testing"
)

func TestTleLoadList(t *testing.T) {
	t.Run("init test", func(t *testing.T) {
		err := loadList("./testdata/tlelist.json")
		if err != nil {
			fmt.Println(err)
			t.Fail()
		}
	})
}
