package gosat

import (
	"errors"
	"math"
)

const (
	twoPi   = math.Pi * 2.0
	deg2rad = math.Pi / 180.0
	x2o3    = 2.0 / 3.0
)

type elsetrec struct {
	error         int
	isInit        bool
	method        rune
	operationmode rune

	/* Near Earth */
	isDeepSpace bool
	aycof       float64
	con41       float64
	cc1         float64
	cc4         float64
	cc5         float64
	d2          float64
	d3          float64
	d4          float64
	delmo       float64
	eta         float64
	argpdot     float64
	omgcof      float64
	sinmao      float64
	t           float64
	t2cof       float64
	t3cof       float64
	t4cof       float64
	t5cof       float64
	x1mth2      float64
	x7thm1      float64
	mdot        float64
	nodedot     float64
	xlcof       float64
	xmcof       float64
	nodecf      float64

	/* Deep Space */
	irez  int
	d2201 float64
	d2211 float64
	d3210 float64
	d3222 float64
	d4410 float64
	d4422 float64
	d5220 float64
	d5232 float64
	d5421 float64
	d5433 float64
	dedt  float64
	del1  float64
	del2  float64
	del3  float64
	didt  float64
	dmdt  float64
	dnodt float64
	domdt float64
	e3    float64
	ee2   float64
	peo   float64
	pgho  float64
	pho   float64
	pinco float64
	plo   float64
	se2   float64
	se3   float64
	sgh2  float64
	sgh3  float64
	sgh4  float64
	sh2   float64
	sh3   float64
	si2   float64
	si3   float64
	sl2   float64
	sl3   float64
	sl4   float64
	gsto  float64
	xfact float64
	xgh2  float64
	xgh3  float64
	xgh4  float64
	xh2   float64
	xh3   float64
	xi2   float64
	xi3   float64
	xl2   float64
	xl3   float64
	xl4   float64
	xlamo float64
	zmol  float64
	zmos  float64
	atime float64
	xli   float64
	xni   float64
	/*  */
	a         float64
	altp      float64
	alta      float64
	rcse      float64
	noUnkozai float64

	gravityVars
	meanVars
	Tle

	// TODO: remove unused fields below
	// Additional elements to capture relevant TLE and object information:
	dia_mm      int     //RSO dia in mm // TODO: long originally
	period_sec  float64 //Period in seconds
	active      bool    //"Active S/C" flag (0=n, 1=y)
	not_orbital bool    //"Orbiting S/C" flag (0=n, 1=y)
	rcs_m2      float64 //"RCS (m^2)" storage
}

type initlVars struct {
	ainv   float64
	ao     float64
	con42  float64
	cosio  float64
	cosio2 float64
	eccsq  float64
	omeosq float64
	posq   float64
	rp     float64
	rteosq float64
	sinio  float64
}

type gravityVars struct {
	tumin         float64
	mu            float64
	radiusearthkm float64
	xke           float64
	j2            float64
	j3            float64
	j4            float64
	j3oj2         float64
}

type meanVars struct {
	am,
	argpm,
	em,
	inclm,
	mm,
	nm,
	nodem float64
}

type commonVars struct {
	cnodm  float64
	cosim  float64
	cosomm float64
	emsq   float64
	gam    float64
	rtemsq float64
	s1     float64
	s2     float64
	s3     float64
	s4     float64
	s5     float64
	s6     float64
	s7     float64
	sinim  float64
	sinomm float64
	snodm  float64
	ss1    float64
	ss2    float64
	ss3    float64
	ss4    float64
	ss5    float64
	ss6    float64
	ss7    float64
	sz1    float64
	sz11   float64
	sz12   float64
	sz13   float64
	sz2    float64
	sz21   float64
	sz22   float64
	sz23   float64
	sz3    float64
	sz31   float64
	sz32   float64
	sz33   float64
	z1     float64
	z11    float64
	z12    float64
	z13    float64
	z2     float64
	z21    float64
	z22    float64
	z23    float64
	z3     float64
	z31    float64
	z32    float64
	z33    float64
}

type dpperVars struct {
	isInit                      bool    // read only
	inclo                       float64 // read only
	ep, inclp, nodep, argpp, mp float64
}

/* -----------------------------------------------------------------------------
*
*                           procedure dpper
*
*  this procedure provides deep space long period periodic contributions
*    to the mean elements.  by design, these periodics are zero at epoch.
*    this used to be dscom which included initialization, but it's really a
*    recurring function.
*
*  author        : david vallado                  719-573-2600   28 jun 2005
*
*  references    :
*    hoots, roehrich, norad spacetrack report #3 1980
*    hoots, norad spacetrack report #6 1986
*    hoots, schumacher and glover 2004
*    vallado, crawford, hujsak, kelso  2006
----------------------------------------------------------------------------*/
func (s *elsetrec) dpper(vars *dpperVars) {
	const (
		zns = 1.19459e-5
		zes = 0.01675
		znl = 1.5835218e-4
		zel = 0.05490
	)

	/* --------------- calculate time varying periodics ----------- */
	zm := s.zmos
	if !vars.isInit {
		zm += zns * s.t
	}

	zf := zm + 2.0*zes*math.Sin(zm)
	sinzf := math.Sin(zf)
	f2 := 0.5*sinzf*sinzf - 0.25
	f3 := -0.5 * sinzf * math.Cos(zf)
	ses := s.se2*f2 + s.se3*f3
	sis := s.si2*f2 + s.si3*f3
	sls := s.sl2*f2 + s.sl3*f3 + s.sl4*sinzf
	sghs := s.sgh2*f2 + s.sgh3*f3 + s.sgh4*sinzf
	shs := s.sh2*f2 + s.sh3*f3
	zm = s.zmol
	if !vars.isInit {
		zm += znl * s.t
	}

	zf = zm + 2.0*zel*math.Sin(zm)
	sinzf = math.Sin(zf)
	f2 = 0.5*sinzf*sinzf - 0.25
	f3 = -0.5 * sinzf * math.Cos(zf)
	sel := s.ee2*f2 + s.e3*f3
	sil := s.xi2*f2 + s.xi3*f3
	sll := s.xl2*f2 + s.xl3*f3 + s.xl4*sinzf
	sghl := s.xgh2*f2 + s.xgh3*f3 + s.xgh4*sinzf
	shll := s.xh2*f2 + s.xh3*f3
	pe := ses + sel
	pinc := sis + sil
	pl := sls + sll
	pgh := sghs + sghl
	ph := shs + shll

	if !vars.isInit {
		pe -= s.peo
		pinc -= s.pinco
		pl -= s.plo
		pgh -= s.pgho
		ph -= s.pho
		vars.inclp += pinc
		vars.ep += pe
		sinip := math.Sin(vars.inclp)
		cosip := math.Cos(vars.inclp)

		/* ----------------- apply periodics directly ------------ */
		//  sgp4fix for lyddane choice
		//  strn3 used original inclination - this is technically feasible
		//  gsfc used perturbed inclination - also technically feasible
		//  probably best to readjust the 0.2 limit value and limit discontinuity
		//  0.2 rad = 11.45916 deg
		//  use next line for original strn3 approach and original inclination
		//  if (inclo >= 0.2)
		//  use next line for gsfc version and perturbed inclination
		if vars.inclp >= 0.2 {
			ph /= sinip
			pgh -= cosip * ph
			vars.argpp += pgh
			vars.nodep += ph
			vars.mp += pl
		} else {
			/* ---- apply periodics with lyddane modification ---- */
			sinop := math.Sin(vars.nodep)
			cosop := math.Cos(vars.nodep)
			alfdp := sinip * sinop
			betdp := sinip * cosop
			dalf := ph*cosop + pinc*cosip*sinop
			dbet := -ph*sinop + pinc*cosip*cosop
			alfdp += dalf
			betdp += dbet
			vars.nodep = math.Mod(vars.nodep, twoPi)
			//  sgp4fix for afspc written intrinsic functions
			// nodep used without a trigonometric function ahead
			if (vars.nodep < 0.0) && (s.operationmode == 'a') {
				vars.nodep += twoPi
			}
			xls := vars.mp + vars.argpp + cosip*vars.nodep
			dls := pl + pgh - pinc*vars.nodep*sinip
			xls += dls
			xnoh := vars.nodep
			vars.nodep = math.Atan2(alfdp, betdp)
			//  sgp4fix for afspc written intrinsic functions
			// nodep used without a trigonometric function ahead
			if (vars.nodep < 0.0) && (s.operationmode == 'a') {
				vars.nodep += twoPi
			}
			if math.Abs(xnoh-vars.nodep) > math.Pi {
				if vars.nodep < xnoh {
					vars.nodep += twoPi
				} else {
					vars.nodep -= twoPi
				}
			}
			vars.mp += pl
			vars.argpp = xls - vars.mp - cosip*vars.nodep
		}
	}
}

/*-----------------------------------------------------------------------------
*
*                           procedure initl
*
*  this procedure initializes the spg4 propagator. all the initialization is
*    consolidated here instead of having multiple loops inside other routines.
*
*  author        : david vallado                  719-573-2600   28 jun 2005
*
*  references    :
*    hoots, roehrich, norad spacetrack report #3 1980
*    hoots, norad spacetrack report #6 1986
*    hoots, schumacher and glover 2004
*    vallado, crawford, hujsak, kelso  2006
----------------------------------------------------------------------------*/
func (s *elsetrec) initl(vars *initlVars) {
	/* ------------- calculate auxillary epoch quantities ---------- */
	vars.eccsq = s.ecco * s.ecco
	vars.omeosq = 1.0 - vars.eccsq
	vars.rteosq = math.Sqrt(vars.omeosq)
	vars.cosio = math.Cos(s.inclo)
	vars.cosio2 = vars.cosio * vars.cosio

	/* ------------------ un-kozai the mean motion ----------------- */
	ak := math.Pow(s.xke/s.no, x2o3)
	d1 := 0.75 * s.j2 * (3.0*vars.cosio2 - 1.0) / (vars.rteosq * vars.omeosq)
	del := d1 / (ak * ak)
	adel := ak * (1.0 - del*del - del*(1.0/3.0+134.0*del*del/81.0))
	del = d1 / (adel * adel)
	s.noUnkozai = s.no / (1.0 + del)

	vars.ao = math.Pow(s.xke/s.noUnkozai, x2o3)
	vars.sinio = math.Sin(s.inclo)
	po := vars.ao * vars.omeosq
	vars.con42 = 1.0 - 5.0*vars.cosio2
	s.con41 = -vars.con42 - vars.cosio2 - vars.cosio2
	vars.ainv = 1.0 / vars.ao
	vars.posq = po * po
	vars.rp = vars.ao * (1.0 - s.ecco)
	s.method = 'n'

	ts70 := s.epoch - 7305.0
	ds70 := math.Floor(ts70 + 1.0e-8)
	tfrac := ts70 - ds70
	// find greenwich location at epoch
	c1 := 1.72027916940703639e-2
	thgr70 := 1.7321343856509374
	fk5r := 5.07551419432269442e-15
	c1p2p := c1 + twoPi
	gsto1 := math.Mod(thgr70+c1*ds70+c1p2p*tfrac+ts70*ts70*fk5r, twoPi)
	if gsto1 < 0.0 {
		gsto1 = gsto1 + twoPi
	}
	s.gsto = gstime(s.epoch + 2433281.5)
}

/*
*
* gstime - find greenwich sidereal time from the julian date
*
 */
func gstime(jdut1 float64) float64 {

	tut1 := (jdut1 - 2451545.0) / 36525.0
	temp := -6.2e-6*tut1*tut1*tut1 +
		0.093104*tut1*tut1 +
		(876600.0*3600+8640184.812866)*tut1 + 67310.54841 // seconds
	temp = math.Mod(temp*deg2rad/240.0, twoPi) //360/86400 = 1/240, to deg, to rad

	// ------------------------ check quadrants ---------------------
	if temp < 0.0 {
		temp += twoPi
	}

	return temp
}

/* -----------------------------------------------------------------------------
*
*                           function getgravconst
*
*  this function gets constants for the propagator. note that mu is identified to
*    facilitiate comparisons with newer models. the common useage is wgs72.
*
*  author        : david vallado                  719-573-2600   21 jul 2006
*
*  references    :
*    norad spacetrack report #3
*    vallado, crawford, hujsak, kelso  2006
--------------------------------------------------------------------------- */
func (s *gravityVars) setGravityVars(whichconst string) {
	switch whichconst {
	// -- wgs-72 low precision str#3 constants --
	case "wgs72old":
		s.mu = 398600.79964        // in km3 / s2
		s.radiusearthkm = 6378.135 // km
		s.xke = 0.0743669161       // reciprocal of tumin
		s.tumin = 1.0 / s.xke
		s.j2 = 0.001082616
		s.j3 = -0.00000253881
		s.j4 = -0.00000165597
		s.j3oj2 = s.j3 / s.j2
	// ------------ wgs-72 constants ------------
	case "wgs72":
		s.mu = 398600.8            // in km3 / s2
		s.radiusearthkm = 6378.135 // km
		s.xke = 60.0 / math.Sqrt(s.radiusearthkm*s.radiusearthkm*s.radiusearthkm/s.mu)
		s.tumin = 1.0 / s.xke
		s.j2 = 0.001082616
		s.j3 = -0.00000253881
		s.j4 = -0.00000165597
		s.j3oj2 = s.j3 / s.j2
	case "wgs84":
		fallthrough
	default:
		// ------------ wgs-84 constants ------------
		s.mu = 398600.5            // in km3 / s2
		s.radiusearthkm = 6378.137 // km
		s.xke = 60.0 / math.Sqrt(s.radiusearthkm*s.radiusearthkm*s.radiusearthkm/s.mu)
		s.tumin = 1.0 / s.xke
		s.j2 = 0.00108262998905
		s.j3 = -0.00000253215306
		s.j4 = -0.00000161098761
		s.j3oj2 = s.j3 / s.j2
	}
}

/*-----------------------------------------------------------------------------
*
*                           procedure dspace
*
*  this procedure provides deep space contributions to mean elements for
*    perturbing third body.  these effects have been averaged over one
*    revolution of the sun and moon.  for earth resonance effects, the
*    effects have been averaged over no revolutions of the satellite.
*    (mean motion)
*
*  author        : david vallado                  719-573-2600   28 jun 2005
*
*  references    :
*    hoots, roehrich, norad spacetrack report #3 1980
*    hoots, norad spacetrack report #6 1986
*    hoots, schumacher and glover 2004
*    vallado, crawford, hujsak, kelso  2006
----------------------------------------------------------------------------*/
func (s *elsetrec) dspace() {
	const (
		fasx2 = 0.13130908
		fasx4 = 2.8843198
		fasx6 = 0.37448087
		g22   = 5.7686396
		g32   = 0.95240898
		g44   = 1.8014998
		g52   = 1.0508330
		g54   = 4.4108898
		rptim = 4.37526908801129966e-3 // this equates to 7.29211514668855e-5 rad/sec
		stepp = 720.0
		stepn = -720.0
		step2 = 259200.0
	)

	/* ----------- calculate deep space resonance effects ----------- */
	theta := math.Mod(s.gsto+s.t*rptim, twoPi)
	s.em += s.dedt * s.t

	s.inclm += s.didt * s.t
	s.argpm += s.domdt * s.t
	s.nodem += s.dnodt * s.t
	s.mm += s.dmdt * s.t

	//   sgp4fix for negative inclinations
	//   the following if statement should be commented out
	//  if (inclm < 0.0)
	// {
	//    inclm = -inclm;
	//    argpm = argpm - pi;
	//    nodem = nodem + pi;
	//  }

	/* - update resonances : numerical (euler-maclaurin) integration - */
	/* ------------------------- epoch restart ----------------------  */
	//   sgp4fix for propagator problems
	//   the following integration works for negative time steps and periods
	//   the specific changes are unknown because the original code was so convoluted

	// sgp4fix take out atime = 0.0 and fix for faster operation
	ft := 0.0
	if s.irez != 0 {
		// sgp4fix streamline check
		if (s.atime == 0.0) || (s.t*s.atime <= 0.0) || (math.Abs(s.t) < math.Abs(s.atime)) {
			s.atime = 0.0
			s.xni = s.noUnkozai
			s.xli = s.xlamo
		}
		// sgp4fix move check outside loop
		var delt float64
		if s.t > 0.0 {
			delt = stepp
		} else {
			delt = stepn
		}

		iretn := 381 // added for do loop
		xndt := 0.0
		xldot := 0.0
		xnddt := 0.0
		xomi := 0.0
		x2li := 0.0
		x2omi := 0.0
		for iretn == 381 {
			/* ------------------- dot terms calculated ------------- */
			/* ----------- near - synchronous resonance terms ------- */
			if s.irez != 2 {
				xndt = s.del1*math.Sin(s.xli-fasx2) + s.del2*math.Sin(2.0*(s.xli-fasx4)) +
					s.del3*math.Sin(3.0*(s.xli-fasx6))
				xldot = s.xni + s.xfact
				xnddt = s.del1*math.Cos(s.xli-fasx2) + 2.0*s.del2*math.Cos(2.0*(s.xli-fasx4)) +
					3.0*s.del3*math.Cos(3.0*(s.xli-fasx6))
				xnddt = xnddt * xldot
			} else {
				/* --------- near - half-day resonance terms -------- */
				xomi = s.argpo + s.argpdot*s.atime
				x2omi = xomi + xomi
				x2li = s.xli + s.xli
				xndt = s.d2201*math.Sin(x2omi+s.xli-g22) + s.d2211*math.Sin(s.xli-g22) +
					s.d3210*math.Sin(xomi+s.xli-g32) + s.d3222*math.Sin(-xomi+s.xli-g32) +
					s.d4410*math.Sin(x2omi+x2li-g44) + s.d4422*math.Sin(x2li-g44) +
					s.d5220*math.Sin(xomi+s.xli-g52) + s.d5232*math.Sin(-xomi+s.xli-g52) +
					s.d5421*math.Sin(xomi+x2li-g54) + s.d5433*math.Sin(-xomi+x2li-g54)
				xldot = s.xni + s.xfact
				xnddt = s.d2201*math.Cos(x2omi+s.xli-g22) + s.d2211*math.Cos(s.xli-g22) +
					s.d3210*math.Cos(xomi+s.xli-g32) + s.d3222*math.Cos(-xomi+s.xli-g32) +
					s.d5220*math.Cos(xomi+s.xli-g52) + s.d5232*math.Cos(-xomi+s.xli-g52) +
					2.0*(s.d4410*math.Cos(x2omi+x2li-g44)+s.d4422*math.Cos(x2li-g44)+
						s.d5421*math.Cos(xomi+x2li-g54)+s.d5433*math.Cos(-xomi+x2li-g54))
				xnddt = xnddt * xldot
			}

			/* ----------------------- integrator ------------------- */
			// sgp4fix move end checks to end of routine
			if math.Abs(s.t-s.atime) >= stepp {
				iretn = 381
			} else { // exit here
				ft = s.t - s.atime
				iretn = 0
			}

			if iretn == 381 {
				s.xli += xldot*delt + xndt*step2
				s.xni += xndt*delt + xnddt*step2
				s.atime += delt
			}
		}

		s.nm = s.xni + xndt*ft + xnddt*ft*ft*0.5
		xl := s.xli + xldot*ft + xndt*ft*ft*0.5
		if s.irez != 1 {
			s.mm = xl - 2.0*s.nodem + 2.0*theta
		} else {
			s.mm = xl - s.nodem - s.argpm + theta
		}
		dndt := s.nm - s.noUnkozai
		s.nm = s.noUnkozai + dndt
	}
}

/*-----------------------------------------------------------------------------
*
*                           procedure dscom
*
*  this procedure provides deep space common items used by both the secular
*    and periodics subroutines.  input is provided as shown. this routine
*    used to be called dpper, but the functions inside weren't well organized.
*
*  author        : david vallado                  719-573-2600   28 jun 2005
*
*  references    :
*    hoots, roehrich, norad spacetrack report #3 1980
*    hoots, norad spacetrack report #6 1986
*    hoots, schumacher and glover 2004
*    vallado, crawford, hujsak, kelso  2006
----------------------------------------------------------------------------*/
func (s *elsetrec) dscom(mv *meanVars, ds *commonVars, tc float64) {
	const (
		zes    = 0.01675
		zel    = 0.05490
		c1ss   = 2.9864797e-6
		c1l    = 4.7968065e-7
		zsinis = 0.39785416
		zcosis = 0.91744867
		zcosgs = 0.1945905
		zsings = -0.98088458
	)

	mv.inclm = s.inclo
	mv.nm = s.noUnkozai
	mv.em = s.ecco
	ds.snodm = math.Sin(s.nodeo)
	ds.cnodm = math.Cos(s.nodeo)
	ds.sinomm = math.Sin(s.argpo)
	ds.cosomm = math.Cos(s.argpo)
	ds.sinim = math.Sin(s.inclo)
	ds.cosim = math.Cos(s.inclo)
	ds.emsq = mv.em * mv.em
	betasq := 1.0 - ds.emsq
	ds.rtemsq = math.Sqrt(betasq)

	/* ----------------- initialize lunar solar terms --------------- */
	s.peo = 0.0
	s.pinco = 0.0
	s.plo = 0.0
	s.pgho = 0.0
	s.pho = 0.0
	day := s.epoch + 18261.5 + tc/1440.0
	xnodce := math.Mod(4.5236020-9.2422029e-4*day, twoPi)
	stem := math.Sin(xnodce)
	ctem := math.Cos(xnodce)
	zcosil := 0.91375164 - 0.03568096*ctem
	zsinil := math.Sqrt(1.0 - zcosil*zcosil)
	zsinhl := 0.089683511 * stem / zsinil
	zcoshl := math.Sqrt(1.0 - zsinhl*zsinhl)
	ds.gam = 5.8351514 + 0.0019443680*day
	zx := 0.39785416 * stem / zsinil
	zy := zcoshl*ctem + 0.91744867*zsinhl*stem
	zx = math.Atan2(zx, zy)
	zx += ds.gam - xnodce
	zcosgl := math.Cos(zx)
	zsingl := math.Sin(zx)

	/* ------------------------- do solar terms --------------------- */
	zcosg := zcosgs
	zsing := zsings
	zcosi := zcosis
	zsini := zsinis
	zcosh := ds.cnodm
	zsinh := ds.snodm
	cc := c1ss
	xnoi := 1.0 / mv.nm

	for lsflg := 1; lsflg <= 2; lsflg++ {
		a1 := zcosg*zcosh + zsing*zcosi*zsinh
		a3 := -zsing*zcosh + zcosg*zcosi*zsinh
		a7 := -zcosg*zsinh + zsing*zcosi*zcosh
		a8 := zsing * zsini
		a9 := zsing*zsinh + zcosg*zcosi*zcosh
		a10 := zcosg * zsini
		a2 := ds.cosim*a7 + ds.sinim*a8
		a4 := ds.cosim*a9 + ds.sinim*a10
		a5 := -ds.sinim*a7 + ds.cosim*a8
		a6 := -ds.sinim*a9 + ds.cosim*a10

		x1 := a1*ds.cosomm + a2*ds.sinomm
		x2 := a3*ds.cosomm + a4*ds.sinomm
		x3 := -a1*ds.sinomm + a2*ds.cosomm
		x4 := -a3*ds.sinomm + a4*ds.cosomm
		x5 := a5 * ds.sinomm
		x6 := a6 * ds.sinomm
		x7 := a5 * ds.cosomm
		x8 := a6 * ds.cosomm

		ds.z31 = 12.0*x1*x1 - 3.0*x3*x3
		ds.z32 = 24.0*x1*x2 - 6.0*x3*x4
		ds.z33 = 12.0*x2*x2 - 3.0*x4*x4
		ds.z1 = 3.0*(a1*a1+a2*a2) + ds.z31*ds.emsq
		ds.z2 = 6.0*(a1*a3+a2*a4) + ds.z32*ds.emsq
		ds.z3 = 3.0*(a3*a3+a4*a4) + ds.z33*ds.emsq
		ds.z11 = -6.0*a1*a5 + ds.emsq*(-24.0*x1*x7-6.0*x3*x5)
		ds.z12 = -6.0*(a1*a6+a3*a5) + ds.emsq*(-24.0*(x2*x7+x1*x8)-6.0*(x3*x6+x4*x5))
		ds.z13 = -6.0*a3*a6 + ds.emsq*(-24.0*x2*x8-6.0*x4*x6)
		ds.z21 = 6.0*a2*a5 + ds.emsq*(24.0*x1*x5-6.0*x3*x7)
		ds.z22 = 6.0*(a4*a5+a2*a6) + ds.emsq*(24.0*(x2*x5+x1*x6)-6.0*(x4*x7+x3*x8))
		ds.z23 = 6.0*a4*a6 + ds.emsq*(24.0*x2*x6-6.0*x4*x8)
		ds.z1 = ds.z1 + ds.z1 + betasq*ds.z31
		ds.z2 = ds.z2 + ds.z2 + betasq*ds.z32
		ds.z3 = ds.z3 + ds.z3 + betasq*ds.z33
		ds.s3 = cc * xnoi
		ds.s2 = -0.5 * ds.s3 / ds.rtemsq
		ds.s4 = ds.s3 * ds.rtemsq
		ds.s1 = -15.0 * mv.em * ds.s4
		ds.s5 = x1*x3 + x2*x4
		ds.s6 = x2*x3 + x1*x4
		ds.s7 = x2*x4 - x1*x3

		/* ----------------------- do lunar terms ------------------- */
		if lsflg == 1 {
			ds.ss1 = ds.s1
			ds.ss2 = ds.s2
			ds.ss3 = ds.s3
			ds.ss4 = ds.s4
			ds.ss5 = ds.s5
			ds.ss6 = ds.s6
			ds.ss7 = ds.s7
			ds.sz1 = ds.z1
			ds.sz2 = ds.z2
			ds.sz3 = ds.z3
			ds.sz11 = ds.z11
			ds.sz12 = ds.z12
			ds.sz13 = ds.z13
			ds.sz21 = ds.z21
			ds.sz22 = ds.z22
			ds.sz23 = ds.z23
			ds.sz31 = ds.z31
			ds.sz32 = ds.z32
			ds.sz33 = ds.z33
			zcosg = zcosgl
			zsing = zsingl
			zcosi = zcosil
			zsini = zsinil
			zcosh = zcoshl*ds.cnodm + zsinhl*ds.snodm
			zsinh = ds.snodm*zcoshl - ds.cnodm*zsinhl
			cc = c1l
		}
	}

	s.zmol = math.Mod(4.7199672+0.22997150*day-ds.gam, twoPi)
	s.zmos = math.Mod(6.2565837+0.017201977*day, twoPi)

	/* ------------------------ do solar terms ---------------------- */
	s.se2 = 2.0 * ds.ss1 * ds.ss6
	s.se3 = 2.0 * ds.ss1 * ds.ss7
	s.si2 = 2.0 * ds.ss2 * ds.sz12
	s.si3 = 2.0 * ds.ss2 * (ds.sz13 - ds.sz11)
	s.sl2 = -2.0 * ds.ss3 * ds.sz2
	s.sl3 = -2.0 * ds.ss3 * (ds.sz3 - ds.sz1)
	s.sl4 = -2.0 * ds.ss3 * (-21.0 - 9.0*ds.emsq) * zes
	s.sgh2 = 2.0 * ds.ss4 * ds.sz32
	s.sgh3 = 2.0 * ds.ss4 * (ds.sz33 - ds.sz31)
	s.sgh4 = -18.0 * ds.ss4 * zes
	s.sh2 = -2.0 * ds.ss2 * ds.sz22
	s.sh3 = -2.0 * ds.ss2 * (ds.sz23 - ds.sz21)

	/* ------------------------ do lunar terms ---------------------- */
	s.ee2 = 2.0 * ds.s1 * ds.s6
	s.e3 = 2.0 * ds.s1 * ds.s7
	s.xi2 = 2.0 * ds.s2 * ds.z12
	s.xi3 = 2.0 * ds.s2 * (ds.z13 - ds.z11)
	s.xl2 = -2.0 * ds.s3 * ds.z2
	s.xl3 = -2.0 * ds.s3 * (ds.z3 - ds.z1)
	s.xl4 = -2.0 * ds.s3 * (-21.0 - 9.0*ds.emsq) * zel
	s.xgh2 = 2.0 * ds.s4 * ds.z32
	s.xgh3 = 2.0 * ds.s4 * (ds.z33 - ds.z31)
	s.xgh4 = -18.0 * ds.s4 * zel
	s.xh2 = -2.0 * ds.s2 * ds.z22
	s.xh3 = -2.0 * ds.s2 * (ds.z23 - ds.z21)
}

/*-----------------------------------------------------------------------------
*
*                           procedure dsinit
*
*  this procedure provides deep space contributions to mean motion dot due
*    to geopotential resonance with half day and one day orbits.
*
*  author        : david vallado                  719-573-2600   28 jun 2005
*
*  references    :
*    hoots, roehrich, norad spacetrack report #3 1980
*    hoots, norad spacetrack report #6 1986
*    hoots, schumacher and glover 2004
*    vallado, crawford, hujsak, kelso  2006
----------------------------------------------------------------------------*/
func (s *elsetrec) dsinit(mv *meanVars, ds *commonVars, iv *initlVars, tc, xpidot float64) {
	const (
		q22    = 1.7891679e-6
		q31    = 2.1460748e-6
		q33    = 2.2123015e-7
		root22 = 1.7891679e-6
		root44 = 7.3636953e-9
		root54 = 2.1765803e-9
		rptim  = 4.37526908801129966e-3 // this equates to 7.29211514668855e-5 rad/sec
		root32 = 3.7393792e-7
		root52 = 1.1428639e-7
		znl    = 1.5835218e-4
		zns    = 1.19459e-5
	)

	aonv := 0.0

	/* -------------------- deep space initialization ------------ */
	if (mv.nm >= 8.26e-3) && (mv.nm <= 9.24e-3) && (mv.em >= 0.5) {
		s.irez = 2
	} else if (mv.nm < 0.0052359877) && (mv.nm > 0.0034906585) {
		s.irez = 1
	} else {
		s.irez = 0
	}

	/* ------------------------ do solar terms ------------------- */
	ses := ds.ss1 * zns * ds.ss5
	sis := ds.ss2 * zns * (ds.sz11 + ds.sz13)
	sls := -zns * ds.ss3 * (ds.sz1 + ds.sz3 - 14.0 - 6.0*ds.emsq)
	sghs := ds.ss4 * zns * (ds.sz31 + ds.sz33 - 6.0)
	shs := -zns * ds.ss2 * (ds.sz21 + ds.sz23)
	// sgp4fix for 180 deg incl
	if (mv.inclm < 5.2359877e-2) || (mv.inclm > math.Pi-5.2359877e-2) {
		shs = 0.0
	}
	if ds.sinim != 0.0 { // TODO: comparing a floating point number with zero
		shs = shs / ds.sinim
	}
	sgs := sghs - ds.cosim*shs

	/* ------------------------- do lunar terms ------------------ */
	s.dedt = ses + ds.s1*znl*ds.s5
	s.didt = sis + ds.s2*znl*(ds.z11+ds.z13)
	s.dmdt = sls - znl*ds.s3*(ds.z1+ds.z3-14.0-6.0*ds.emsq)
	sghl := ds.s4 * znl * (ds.z31 + ds.z33 - 6.0)
	shll := -znl * ds.s2 * (ds.z21 + ds.z23)
	// sgp4fix for 180 deg incl
	if (mv.inclm < 5.2359877e-2) || (mv.inclm > math.Pi-5.2359877e-2) {
		shll = 0.0
	}
	s.domdt = sgs + sghl
	s.dnodt = shs
	if ds.sinim != 0.0 { // TODO: comparing a floating point number with zero
		s.domdt -= ds.cosim / ds.sinim * shll
		s.dnodt += shll / ds.sinim
	}

	/* ----------- calculate deep space resonance effects -------- */
	theta := math.Mod(s.gsto+tc*rptim, twoPi)
	mv.em += s.dedt * s.t
	mv.inclm += s.didt * s.t // TODO: these calculations aren't used further
	mv.argpm += s.domdt * s.t
	mv.nodem += s.dnodt * s.t
	mv.mm += s.dmdt * s.t
	//   sgp4fix for negative inclinations
	//   the following if statement should be commented out
	// if (inclm < 0.0)
	//  {
	//    inclm  = -inclm;
	//    argpm  = argpm - pi;
	//    nodem = nodem + pi;
	//  }

	/* -------------- initialize the resonance terms ------------- */
	if s.irez != 0 {
		aonv = math.Pow(mv.nm/s.xke, x2o3)

		/* ---------- geopotential resonance for 12 hour orbits ------ */
		if s.irez == 2 {
			cosisq := ds.cosim * ds.cosim
			emo := mv.em
			mv.em = s.ecco
			emsqo := ds.emsq
			ds.emsq = iv.eccsq
			eoc := mv.em * ds.emsq
			g201 := -0.306 - (mv.em-0.64)*0.440

			var g211, g310, g322, g410, g422, g520, g521, g532, g533 float64

			if mv.em <= 0.65 {
				g211 = 3.616 - 13.2470*mv.em + 16.2900*ds.emsq
				g310 = -19.302 + 117.3900*mv.em - 228.4190*ds.emsq + 156.5910*eoc
				g322 = -18.9068 + 109.7927*mv.em - 214.6334*ds.emsq + 146.5816*eoc
				g410 = -41.122 + 242.6940*mv.em - 471.0940*ds.emsq + 313.9530*eoc
				g422 = -146.407 + 841.8800*mv.em - 1629.014*ds.emsq + 1083.4350*eoc
				g520 = -532.114 + 3017.977*mv.em - 5740.032*ds.emsq + 3708.2760*eoc
			} else {
				g211 = -72.099 + 331.819*mv.em - 508.738*ds.emsq + 266.724*eoc
				g310 = -346.844 + 1582.851*mv.em - 2415.925*ds.emsq + 1246.113*eoc
				g322 = -342.585 + 1554.908*mv.em - 2366.899*ds.emsq + 1215.972*eoc
				g410 = -1052.797 + 4758.686*mv.em - 7193.992*ds.emsq + 3651.957*eoc
				g422 = -3581.690 + 16178.110*mv.em - 24462.770*ds.emsq + 12422.520*eoc
				if mv.em > 0.715 {
					g520 = -5149.66 + 29936.92*mv.em - 54087.36*ds.emsq + 31324.56*eoc
				} else {
					g520 = 1464.74 - 4664.75*mv.em + 3763.64*ds.emsq
				}
			}
			if mv.em < 0.7 {
				g533 = -919.22770 + 4988.6100*mv.em - 9064.7700*ds.emsq + 5542.21*eoc
				g521 = -822.71072 + 4568.6173*mv.em - 8491.4146*ds.emsq + 5337.524*eoc
				g532 = -853.66600 + 4690.2500*mv.em - 8624.7700*ds.emsq + 5341.4*eoc
			} else {
				g533 = -37995.780 + 161616.52*mv.em - 229838.20*ds.emsq + 109377.94*eoc
				g521 = -51752.104 + 218913.95*mv.em - 309468.16*ds.emsq + 146349.42*eoc
				g532 = -40023.880 + 170470.89*mv.em - 242699.48*ds.emsq + 115605.82*eoc
			}

			sini2 := ds.sinim * ds.sinim
			f220 := 0.75 * (1.0 + 2.0*ds.cosim + cosisq)
			f221 := 1.5 * sini2
			f321 := 1.875 * ds.sinim * (1.0 - 2.0*ds.cosim - 3.0*cosisq)
			f322 := -1.875 * ds.sinim * (1.0 + 2.0*ds.cosim - 3.0*cosisq)
			f441 := 35.0 * sini2 * f220
			f442 := 39.3750 * sini2 * sini2
			f522 := 9.84375 * ds.sinim * (sini2*(1.0-2.0*ds.cosim-5.0*cosisq) +
				0.33333333*(-2.0+4.0*ds.cosim+6.0*cosisq))
			f523 := ds.sinim * (4.92187512*sini2*(-2.0-4.0*ds.cosim+10.0*cosisq) +
				6.56250012*(1.0+2.0*ds.cosim-3.0*cosisq))
			f542 := 29.53125 * ds.sinim *
				(2.0 - 8.0*ds.cosim + cosisq*(-12.0+8.0*ds.cosim+10.0*cosisq))
			f543 := 29.53125 * ds.sinim *
				(-2.0 - 8.0*ds.cosim + cosisq*(12.0+8.0*ds.cosim-10.0*cosisq))
			xno2 := mv.nm * mv.nm
			ainv2 := aonv * aonv
			temp1 := 3.0 * xno2 * ainv2
			temp := temp1 * root22
			s.d2201 = temp * f220 * g201
			s.d2211 = temp * f221 * g211
			temp1 = temp1 * aonv
			temp = temp1 * root32
			s.d3210 = temp * f321 * g310
			s.d3222 = temp * f322 * g322
			temp1 = temp1 * aonv
			temp = 2.0 * temp1 * root44
			s.d4410 = temp * f441 * g410
			s.d4422 = temp * f442 * g422
			temp1 = temp1 * aonv
			temp = temp1 * root52
			s.d5220 = temp * f522 * g520
			s.d5232 = temp * f523 * g532
			temp = 2.0 * temp1 * root54
			s.d5421 = temp * f542 * g521
			s.d5433 = temp * f543 * g533
			s.xlamo = math.Mod(s.mo+s.nodeo+s.nodeo-theta-theta, twoPi)
			s.xfact = s.mdot + s.dmdt + 2.0*(s.nodedot+s.dnodt-rptim) - s.noUnkozai
			mv.em = emo
			ds.emsq = emsqo
		}

		/* ---------------- synchronous resonance terms -------------- */
		if s.irez == 1 {
			g200 := 1.0 + ds.emsq*(-2.5+0.8125*ds.emsq)
			g310 := 1.0 + 2.0*ds.emsq
			g300 := 1.0 + ds.emsq*(-6.0+6.60937*ds.emsq)
			f220 := 0.75 * (1.0 + ds.cosim) * (1.0 + ds.cosim)
			f311 := 0.9375*ds.sinim*ds.sinim*(1.0+3.0*ds.cosim) - 0.75*(1.0+ds.cosim)
			f330 := 1.0 + ds.cosim
			f330 = 1.875 * f330 * f330 * f330
			s.del1 = 3.0 * mv.nm * mv.nm * aonv * aonv
			s.del2 = 2.0 * s.del1 * f220 * g200 * q22
			s.del3 = 3.0 * s.del1 * f330 * g300 * q33 * aonv
			s.del1 = s.del1 * f311 * g310 * q31 * aonv
			s.xlamo = math.Mod(s.mo+s.nodeo+s.argpo-theta, twoPi)
			s.xfact = s.mdot + xpidot - rptim + s.dmdt + s.domdt + s.dnodt - s.noUnkozai
		}

		/* ------------ for sgp4, initialize the integrator ---------- */
		s.xli = s.xlamo
		s.xni = s.noUnkozai
		s.atime = 0.0
		mv.nm = s.noUnkozai
	}
}

func (s *elsetrec) sgp4init(whichconst string, t *Tle, opsmode rune) error {
	const temp4 = 1.5e-12

	/* ----------- set all near earth variables to zero ------------ */
	s.isDeepSpace = false
	s.method = 'n'
	s.aycof = 0.0
	s.con41 = 0.0
	s.cc1 = 0.0
	s.cc4 = 0.0
	s.cc5 = 0.0
	s.d2 = 0.0
	s.d3 = 0.0
	s.d4 = 0.0
	s.delmo = 0.0
	s.eta = 0.0
	s.argpdot = 0.0
	s.omgcof = 0.0
	s.sinmao = 0.0
	s.t = 0.0
	s.t2cof = 0.0
	s.t3cof = 0.0
	s.t4cof = 0.0
	s.t5cof = 0.0
	s.x1mth2 = 0.0
	s.x7thm1 = 0.0
	s.mdot = 0.0
	s.nodedot = 0.0
	s.xlcof = 0.0
	s.xmcof = 0.0
	s.nodecf = 0.0

	/* ----------- set all deep space variables to zero ------------ */
	s.irez = 0
	s.d2201 = 0.0
	s.d2211 = 0.0
	s.d3210 = 0.0
	s.d3222 = 0.0
	s.d4410 = 0.0
	s.d4422 = 0.0
	s.d5220 = 0.0
	s.d5232 = 0.0
	s.d5421 = 0.0
	s.d5433 = 0.0
	s.dedt = 0.0
	s.del1 = 0.0
	s.del2 = 0.0
	s.del3 = 0.0
	s.didt = 0.0
	s.dmdt = 0.0
	s.dnodt = 0.0
	s.domdt = 0.0
	s.e3 = 0.0
	s.ee2 = 0.0
	s.peo = 0.0
	s.pgho = 0.0
	s.pho = 0.0
	s.pinco = 0.0
	s.plo = 0.0
	s.se2 = 0.0
	s.se3 = 0.0
	s.sgh2 = 0.0
	s.sgh3 = 0.0
	s.sgh4 = 0.0
	s.sh2 = 0.0
	s.sh3 = 0.0
	s.si2 = 0.0
	s.si3 = 0.0
	s.sl2 = 0.0
	s.sl3 = 0.0
	s.sl4 = 0.0
	s.gsto = 0.0
	s.xfact = 0.0
	s.xgh2 = 0.0
	s.xgh3 = 0.0
	s.xgh4 = 0.0
	s.xh2 = 0.0
	s.xh3 = 0.0
	s.xi2 = 0.0
	s.xi3 = 0.0
	s.xl2 = 0.0
	s.xl3 = 0.0
	s.xl4 = 0.0
	s.xlamo = 0.0
	s.zmol = 0.0
	s.zmos = 0.0
	s.atime = 0.0
	s.xli = 0.0
	s.xni = 0.0

	/* ------------------------ earth constants ----------------------- */
	s.setGravityVars(whichconst)
	s.copy(t)

	s.error = 0
	s.operationmode = opsmode

	s.am = 0.0
	s.em = 0.0
	s.inclm = 0.0
	s.nodem = 0.0
	s.argpm = 0.0
	s.mm = 0.0
	s.nm = 0.0

	ss := 78.0/s.radiusearthkm + 1.0
	qzms2t := (120.0 - 78.0) / s.radiusearthkm
	qzms2t *= qzms2t
	qzms2t *= qzms2t // pow 4

	s.isInit = true
	s.t = 0.0

	var iv initlVars
	s.initl(&iv)

	s.a = math.Pow(s.noUnkozai*s.tumin, -x2o3)
	s.alta = s.a*(1.0+s.ecco) - 1.0
	s.altp = s.a*(1.0-s.ecco) - 1.0
	s.error = 0

	if (iv.omeosq >= 0.0) || (s.noUnkozai >= 0.0) {
		if iv.rp < (220.0/s.radiusearthkm + 1.0) {
			s.isDeepSpace = true
		} else {
			s.isDeepSpace = false
		}
		sfour := ss
		qzms24 := qzms2t
		perige := (iv.rp - 1.0) * s.radiusearthkm

		/* - for perigees below 156 km, s and qoms2t are altered - */
		if perige < 156.0 {
			sfour = perige - 78.0
			if perige < 98.0 {
				sfour = 20.0
			}
			qzms24temp := (120.0 - sfour) / s.radiusearthkm
			qzms24 = qzms24temp * qzms24temp * qzms24temp * qzms24temp
			sfour = sfour/s.radiusearthkm + 1.0
		}
		pinvsq := 1.0 / iv.posq
		tsi := 1.0 / (iv.ao - sfour)
		s.eta = iv.ao * s.ecco * tsi
		etasq := s.eta * s.eta
		eeta := s.ecco * s.eta
		psisq := math.Abs(1.0 - etasq)
		tsisq := tsi * tsi
		coef := qzms24 * tsisq * tsisq
		coef1 := coef / math.Pow(psisq, 3.5)
		cc2 := coef1 * s.noUnkozai * (iv.ao*(1.0+1.5*etasq+eeta*(4.0+etasq)) +
			0.375*s.j2*tsi/psisq*s.con41*(8.0+3.0*etasq*(8.0+etasq)))
		s.cc1 = s.bstar * cc2
		cc3 := 0.0
		if s.ecco > 1.0e-4 {
			cc3 = -2.0 * coef * tsi * s.j3oj2 * s.noUnkozai * iv.sinio / s.ecco
		}
		s.x1mth2 = 1.0 - iv.cosio2
		s.cc4 = 2.0 * s.noUnkozai * coef1 * iv.ao * iv.omeosq * (s.eta*(2.0+0.5*etasq) +
			s.ecco*(0.5+2.0*etasq) -
			s.j2*tsi/(iv.ao*psisq)*(-3.0*s.con41*(1.0-2.0*eeta+etasq*(1.5-0.5*eeta))+
				0.75*s.x1mth2*(2.0*etasq-eeta*(1.0+etasq))*math.Cos(2.0*s.argpo)))
		s.cc5 = 2.0 * coef1 * iv.ao * iv.omeosq * (1.0 + 2.75*(etasq+eeta) + eeta*etasq)
		cosio4 := iv.cosio2 * iv.cosio2
		temp1 := 1.5 * s.j2 * pinvsq * s.noUnkozai
		temp2 := 0.5 * temp1 * s.j2 * pinvsq
		temp3 := -0.46875 * s.j4 * pinvsq * pinvsq * s.noUnkozai
		s.mdot = s.noUnkozai + 0.5*temp1*iv.rteosq*s.con41 +
			0.0625*temp2*iv.rteosq*(13.0-78.0*iv.cosio2+137.0*cosio4)
		s.argpdot = -0.5*temp1*iv.con42 + 0.0625*temp2*(7.0-114.0*iv.cosio2+395.0*cosio4) +
			temp3*(3.0-36.0*iv.cosio2+49.0*cosio4)
		xhdot1 := -temp1 * iv.cosio
		s.nodedot = xhdot1 + (0.5*temp2*(4.0-19.0*iv.cosio2)+2.0*temp3*(3.0-7.0*iv.cosio2))*iv.cosio
		xpidot := s.argpdot + s.nodedot
		s.omgcof = s.bstar * cc3 * math.Cos(s.argpo)
		if s.ecco > 1.0e-4 {
			s.xmcof = -x2o3 * coef * s.bstar / eeta
		} else {
			s.xmcof = 0.0
		}
		s.nodecf = 3.5 * iv.omeosq * xhdot1 * s.cc1
		s.t2cof = 1.5 * s.cc1

		if math.Abs(iv.cosio+1.0) > 1.5e-12 {
			s.xlcof = -0.25 * s.j3oj2 * iv.sinio * (3.0 + 5.0*iv.cosio) / (1.0 + iv.cosio)
		} else {
			s.xlcof = -0.25 * s.j3oj2 * iv.sinio * (3.0 + 5.0*iv.cosio) / temp4
		}
		s.aycof = -0.5 * s.j3oj2 * iv.sinio
		delmotemp := 1.0 + s.eta*math.Cos(s.mo)
		s.delmo = delmotemp * delmotemp * delmotemp
		s.sinmao = math.Sin(s.mo)
		s.x7thm1 = 7.0*iv.cosio2 - 1.0

		/* --------------- deep space initialization ------------- */
		if (twoPi / s.noUnkozai) >= 225.0 {
			s.method = 'd'
			s.isDeepSpace = true
			tc := 0.0 // TODO: remove tc
			var ds commonVars
			var mv meanVars
			s.dscom(&mv, &ds, tc)

			var dpperVars dpperVars
			dpperVars.isInit = s.isInit
			dpperVars.inclo = s.inclo
			dpperVars.ep = s.ecco
			dpperVars.inclp = s.inclo
			dpperVars.nodep = s.nodeo
			dpperVars.argpp = s.argpo
			dpperVars.mp = s.mo

			s.dpper(&dpperVars)
			s.ecco = dpperVars.ep
			s.inclo = dpperVars.inclp
			s.nodeo = dpperVars.nodep
			s.argpo = dpperVars.argpp
			s.mo = dpperVars.mp

			mv.argpm = 0.0
			mv.nodem = 0.0
			mv.mm = 0.0
			s.dsinit(&mv, &ds, &iv, tc, xpidot)
		}

		/* ----------- set variables if not deep space ----------- */
		if !s.isDeepSpace {
			cc1sq := s.cc1 * s.cc1
			s.d2 = 4.0 * iv.ao * tsi * cc1sq
			temp := s.d2 * tsi * s.cc1 / 3.0
			s.d3 = (17.0*iv.ao + sfour) * temp
			s.d4 = 0.5 * temp * iv.ao * tsi * (221.0*iv.ao + 31.0*sfour) * s.cc1
			s.t3cof = s.d2 + 2.0*cc1sq
			s.t4cof = 0.25 * (3.0*s.d3 + s.cc1*(12.0*s.d2+10.0*cc1sq))
			s.t5cof = 0.2 * (3.0*s.d4 + 12.0*s.cc1*s.d3 + 6.0*s.d2*s.d2 + 15.0*cc1sq*(2.0*s.d2+cc1sq))
		}
	}

	/* finally propogate to zero epoch to initialize all others. */
	// sgp4fix take out check to let satellites process until they are actually below earth surface
	//       if(satrec.error == 0)
	var r [6]float64
	err := s.sgp4(0.0, r[0:3], r[3:])

	s.isInit = false

	return err
}

/*-----------------------------------------------------------------------------
*
*                             procedure sgp4
*
*  this procedure is the sgp4 prediction model from space command. this is an
*    updated and combined version of sgp4 and sdp4, which were originally
*    published separately in spacetrack report #3. this version follows the
*    methodology from the aiaa paper (2006) describing the history and
*    development of the code.
*
*  author        : david vallado                  719-573-2600   28 jun 2005
*
*  inputs        :
*    satrec	 - initialised structure from sgp4init() call.
*    tsince	 - time since epoch (minutes)
*
*  outputs       :
*    r           - position vector                     km
*    v           - velocity                            km/sec
*  return code - non-zero on error.
*                   1 - mean elements, ecc >= 1.0 or ecc < -0.001 or a < 0.95 er
*                   2 - mean motion less than 0.0
*                   3 - pert elements, ecc < 0.0  or  ecc > 1.0
*                   4 - semi-latus rectum < 0.0
*                   5 - epoch elements are sub-orbital
*                   6 - satellite has decayed
*
*  references    :
*    hoots, roehrich, norad spacetrack report #3 1980
*    hoots, norad spacetrack report #6 1986
*    hoots, schumacher and glover 2004
*    vallado, crawford, hujsak, kelso  2006
----------------------------------------------------------------------------*/

func (s *elsetrec) sgp4(tsince float64, r []float64, v []float64) error {
	/* ------------------ set mathematical constants --------------- */
	// sgp4fix divisor for divide by zero check on inclination
	// the old check used 1.0 + cos(pi-1.0e-9), but then compared it to
	// 1.5 e-12, so the threshold was changed to 1.5e-12 for consistency
	const temp4 = 1.5e-12

	// var meanVars meanVars // TODO: remove local this var
	vkmpersec := s.radiusearthkm * s.xke / 60.0

	/* --------------------- clear sgp4 error flag ----------------- */
	s.t = tsince
	s.error = 0

	/* ------- update for secular gravity and atmospheric drag ----- */
	xmdf := s.mo + s.mdot*s.t
	argpdf := s.argpo + s.argpdot*s.t
	nodedf := s.nodeo + s.nodedot*s.t
	s.argpm = argpdf
	s.mm = xmdf
	t2 := s.t * s.t
	s.nodem = nodedf + s.nodecf*t2
	tempa := 1.0 - s.cc1*s.t
	tempe := s.bstar * s.cc4 * s.t
	templ := s.t2cof * t2

	if !s.isDeepSpace {
		delomg := s.omgcof * s.t
		// sgp4fix use mutliply for speed instead of pow
		delmtemp := 1.0 + s.eta*math.Cos(xmdf)
		delm := s.xmcof * (delmtemp*delmtemp*delmtemp - s.delmo)
		temp := delomg + delm
		s.mm = xmdf + temp
		s.argpm = argpdf - temp
		t3 := t2 * s.t
		t4 := t3 * s.t
		tempa -= s.d2*t2 + s.d3*t3 + s.d4*t4
		tempe += s.bstar * s.cc5 * (math.Sin(s.mm) - s.sinmao)
		templ += s.t3cof*t3 + t4*(s.t4cof+s.t*s.t5cof)
	}

	s.nm = s.noUnkozai
	s.em = s.ecco
	s.inclm = s.inclo
	if s.method == 'd' {
		s.dspace()
	}

	if s.nm <= 0.0 {
		s.error = 2
		// sgp4fix add return
		return errors.New("nm <= 0.0")
	}
	s.am = math.Pow(s.xke/s.nm, x2o3) * tempa * tempa
	s.nm = s.xke / math.Pow(s.am, 1.5)
	s.em -= tempe

	// fix tolerance for error recognition
	// sgp4fix am is fixed from the previous nm check
	if (s.em >= 1.0) || (s.em < -0.001) /* || (am < 0.95)*/ {
		s.error = 1
		// sgp4fix to return if there is an error in eccentricity
		return errors.New("(em >= 1.0) || (em < -0.001)")
	}
	// sgp4fix fix tolerance to avoid a divide by zero
	if s.em < 1.0e-6 {
		s.em = 1.0e-6
	}
	s.mm += s.noUnkozai * templ
	xlm := s.mm + s.argpm + s.nodem
	emsq := s.em * s.em
	temp := 1.0 - emsq

	s.nodem = math.Mod(s.nodem, twoPi)
	s.argpm = math.Mod(s.argpm, twoPi)
	xlm = math.Mod(xlm, twoPi)
	s.mm = math.Mod(xlm-s.argpm-s.nodem, twoPi)

	// sgp4fix recover singly averaged mean elements
	// s.am = s.am // TODO
	// s.em = s.em
	// s.inclm = s.inclm
	// s.nodem = s.nodem
	// s.argpm = s.argpm
	// s.mm = s.mm
	// s.nm = s.nm

	/* ----------------- compute extra mean quantities ------------- */
	sinim := math.Sin(s.inclm)
	cosim := math.Cos(s.inclm)

	/* -------------------- add lunar-solar periodics -------------- */
	ep := s.em
	xincp := s.inclm
	argpp := s.argpm
	nodep := s.nodem
	mp := s.mm
	sinip := sinim
	cosip := cosim

	var dpperVars dpperVars
	dpperVars.isInit = false
	dpperVars.inclo = s.inclo
	dpperVars.ep = ep
	dpperVars.inclp = xincp
	dpperVars.nodep = nodep
	dpperVars.argpp = argpp
	dpperVars.mp = mp

	if s.method == 'd' {
		s.dpper(&dpperVars)
		ep = dpperVars.ep
		xincp = dpperVars.inclp
		nodep = dpperVars.nodep
		argpp = dpperVars.argpp
		mp = dpperVars.mp

		if xincp < 0.0 {
			xincp = -xincp
			nodep = nodep + math.Pi
			argpp = argpp - math.Pi
		}
		if (ep < 0.0) || (ep > 1.0) {
			s.error = 3
			// sgp4fix add return
			return errors.New("(ep < 0.0) || (ep > 1.0)")
		}
	}

	/* -------------------- long period periodics ------------------ */
	if s.method == 'd' {
		sinip = math.Sin(xincp)
		cosip = math.Cos(xincp)
		s.aycof = -0.5 * s.j3oj2 * sinip
		// sgp4fix for divide by zero for xincp = 180 deg
		if math.Abs(cosip+1.0) > 1.5e-12 {
			s.xlcof = -0.25 * s.j3oj2 * sinip * (3.0 + 5.0*cosip) / (1.0 + cosip)
		} else {
			s.xlcof = -0.25 * s.j3oj2 * sinip * (3.0 + 5.0*cosip) / temp4
		}
	}
	axnl := ep * math.Cos(argpp)
	temp = 1.0 / (s.am * (1.0 - ep*ep))
	aynl := ep*math.Sin(argpp) + temp*s.aycof
	xl := mp + argpp + nodep + temp*s.xlcof*axnl

	/* --------------------- solve kepler's equation --------------- */
	u := math.Mod(xl-nodep, twoPi)
	eo1 := u
	tem5 := 9999.9
	ktr := 1
	//   sgp4fix for kepler iteration
	//   the following iteration needs better limits on corrections
	sineo1 := 0.0
	coseo1 := 0.0
	for (math.Abs(tem5) >= 1.0e-12) && (ktr <= 10) {
		sineo1 = math.Sin(eo1)
		coseo1 = math.Cos(eo1)
		tem5 = 1.0 - coseo1*axnl - sineo1*aynl
		tem5 = (u - aynl*coseo1 + axnl*sineo1 - eo1) / tem5
		if math.Abs(tem5) >= 0.95 {
			if tem5 > 0.0 {
				tem5 = 0.95
			} else {
				tem5 = -0.95
			}
		}
		eo1 = eo1 + tem5
		ktr = ktr + 1
	}

	/* ------------- short period preliminary quantities ----------- */
	ecose := axnl*coseo1 + aynl*sineo1
	esine := axnl*sineo1 - aynl*coseo1
	el2 := axnl*axnl + aynl*aynl
	pl := s.am * (1.0 - el2)
	if pl < 0.0 {
		s.error = 4
		// sgp4fix add return
		return errors.New("pl < 0.0")
	}
	rl := s.am * (1.0 - ecose)
	rdotl := math.Sqrt(s.am) * esine / rl
	rvdotl := math.Sqrt(pl) / rl
	betal := math.Sqrt(1.0 - el2)
	temp = esine / (1.0 + betal)
	sinu := s.am / rl * (sineo1 - aynl - axnl*temp)
	cosu := s.am / rl * (coseo1 - axnl + aynl*temp)
	su := math.Atan2(sinu, cosu)
	sin2u := (cosu + cosu) * sinu
	cos2u := 1.0 - 2.0*sinu*sinu
	temp = 1.0 / pl
	temp1 := 0.5 * s.j2 * temp
	temp2 := temp1 * temp

	/* -------------- update for short period periodics ------------ */
	if s.method == 'd' {
		cosisq := cosip * cosip
		s.con41 = 3.0*cosisq - 1.0
		s.x1mth2 = 1.0 - cosisq
		s.x7thm1 = 7.0*cosisq - 1.0
	}
	mrt := rl*(1.0-1.5*temp2*betal*s.con41) + 0.5*temp1*s.x1mth2*cos2u
	su = su - 0.25*temp2*s.x7thm1*sin2u
	xnode := nodep + 1.5*temp2*cosip*sin2u
	xinc := xincp + 1.5*temp2*cosip*sinip*cos2u
	mvt := rdotl - s.nm*temp1*s.x1mth2*sin2u/s.xke
	rvdot := rvdotl + s.nm*temp1*(s.x1mth2*cos2u+1.5*s.con41)/s.xke

	/* --------------------- orientation vectors ------------------- */
	sinsu := math.Sin(su)
	cossu := math.Cos(su)
	snod := math.Sin(xnode)
	cnod := math.Cos(xnode)
	sini := math.Sin(xinc)
	cosi := math.Cos(xinc)
	xmx := -snod * cosi
	xmy := cnod * cosi
	ux := xmx*sinsu + cnod*cossu
	uy := xmy*sinsu + snod*cossu
	uz := sini * sinsu
	vx := xmx*cossu - cnod*sinsu
	vy := xmy*cossu - snod*sinsu
	vz := sini * cossu

	/* --------- position and velocity (in km and km/sec) ---------- */
	r[0] = mrt * ux * s.radiusearthkm
	r[1] = mrt * uy * s.radiusearthkm
	r[2] = mrt * uz * s.radiusearthkm
	v[0] = (mvt*ux + rvdot*vx) * vkmpersec
	v[1] = (mvt*uy + rvdot*vy) * vkmpersec
	v[2] = (mvt*uz + rvdot*vz) * vkmpersec

	// sgp4fix for decaying satellites
	if mrt < 1.0 {
		s.error = 6
		return errors.New("mrt < 1.0")
	}

	return nil
}
