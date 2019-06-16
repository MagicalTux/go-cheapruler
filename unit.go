package cheapruler

type Unit float64

const (
	Kilometers    Unit = 1.0
	Kilometres    Unit = 1.0
	Miles         Unit = 1000 / 1609.344
	Nauticalmiles Unit = 1000 / 1852
	Meters        Unit = 1000
	Metres        Unit = 1000
	Yards         Unit = 1000 / 0.9144
	Feet          Unit = 1000 / 0.3048
	Inches        Unit = 1000 / 0.0254
)
