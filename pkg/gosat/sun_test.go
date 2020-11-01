package gosat

import (
	"math"
	"testing"
	"time"
)

const (
	azmTest   = (154.384634974813025110051967203617095947265625) * degToRad
	elvTest   = (5.8779458855349417945035384036600589752197265625) * degToRad
	tolerance = 20.0
)

func TestSunAe(t *testing.T) {
	// 01.01.2018 06:00:00 UTC
	time := time.Unix(1514786400, 0).UTC()
	latitude := 58.0 * degToRad
	longitude := 63.0 * degToRad
	sunLatitude, sunLongitude, _ := SunLlr(time)

	t.Log(time)
	t.Log("Sun location")
	t.Logf("\tLat: %+6.2f\n", sunLatitude*radToDeg)
	t.Logf("\tLon: %+6.2f\n\n", sunLongitude*radToDeg)
	azimuth, elevation := sunAe(time, latitude, longitude, true)

	dazm := math.Abs(azimuth - azmTest)
	delv := math.Abs(elevation - elvTest)

	eazm := math.Max(math.Abs(azimuth), math.Abs(azmTest))
	eazm = (math.Nextafter(eazm, math.MaxFloat64) - eazm) * tolerance
	eelv := math.Max(math.Abs(elevation), math.Abs(elvTest))
	eelv = (math.Nextafter(eelv, math.MaxFloat64) - eelv) * tolerance

	if dazm > eazm || delv > eelv {
		t.Errorf("dazm: %0.3e tolerance: %0.3e\n", dazm, eazm)
		t.Errorf("eazm: %0.3e tolerance: %0.3e\n", delv, eelv)
	}
}
