package gosat

import (
	"encoding/json"
	"fmt"
	"io"
	"io/ioutil"
	"log"
	"net/http"
	"os"
	"path/filepath"
	"strings"
)

type tleList struct {
	TleList []string `json:"tleList"`
}

type writeCounter struct {
	Total uint64
}

func (wc *writeCounter) Write(p []byte) (int, error) {
	n := len(p)
	wc.Total += uint64(n)
	wc.printProgress()
	return n, nil
}

func (wc writeCounter) printProgress() {
	//clear line
	fmt.Printf("\r%s", strings.Repeat(" ", 50))

	var suffix string
	var n uint64
	if wc.Total >= 0x40000000 {
		suffix = "G"
		n = wc.Total / 0x40000000
	} else if wc.Total >= 0x00100000 {
		suffix = "M"
		n = wc.Total / 0x00100000
	} else if wc.Total >= 0x00000400 {
		suffix = "k"
		n = wc.Total / 0x00000400
	} else {
		suffix = ""
		n = wc.Total
	}

	fmt.Printf("\rDownloading %d%sB complete", n, suffix)
}

func loadList(path string) error {
	jsonFile, err := os.Open(path)
	if err != nil {
		fmt.Println(err)
		return err
	}
	defer jsonFile.Close()
	bytes, _ := ioutil.ReadAll(jsonFile)

	var list tleList
	json.Unmarshal(bytes, &list)

	dir, err := os.UserCacheDir()
	if err != nil {
		log.Fatal(err)
	}
	dir = filepath.Join(dir, "gosat")

	for i := 0; i < len(list.TleList); i++ {
		fmt.Println(list.TleList[i])
		file := filepath.Base(list.TleList[i])
		file = filepath.Join(dir, file)
		err := downloadFile(list.TleList[i], file)
		if err != nil {
			log.Fatal(err)
			return err
		}
	}

	return nil
}

func downloadFile(url string, filepath string) error {
	out, err := os.Create(filepath)
	if err != nil {
		return err
	}
	defer out.Close()

	resp, err := http.Get(url)
	if err != nil {
		return err
	}
	defer resp.Body.Close()

	counter := &writeCounter{}
	_, err = io.Copy(out, io.TeeReader(resp.Body, counter))
	if err != nil {
		return err
	}
	fmt.Println()

	return nil
}
