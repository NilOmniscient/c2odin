package c2odin

// Structs  and applicable functions that are used in the lib that aren't shapes.
Simplex :: struct {
  verts: [4]SV,
  div: f32,
  count: int,
}
SV :: struct { // SimplexVert
  sA, sB, p: [2]f32,
  u: f32,
  iA, iB: int
}
Proxy :: struct {
  radius: f32,
  count: int,
  verts: [MAX_POLYGON_VERTS][2]f32,
}
HalfSpace :: struct {
  n: [2]f32,  // Normal
  d: f32,   // Distance
}
Matrix :: [2][2]f32

MakeProxy :: proc(shape: ^Shape) -> Proxy{
  p := Proxy{}
  switch shape.Type {
    case Circle:
      c := cast(^Circle)shape
      p.radius = c.Radius
      p.count = 1
      p.verts[0] = c.Origin
    case AABB:
      bb := cast(^AABB)shape
      p.radius = 0
      p.count = 4
      BBVerts(bb, p.verts[:])
    case Capsule:
      c := cast(^Capsule)shape
      p.radius = c.Radius
      p.count = 2
      p.verts[0] = c.Origin
      p.verts[1] = c.Origin2
    case Polygon:
      poly := cast(^Polygon)shape
      p.radius = 0
      p.count = poly.Count
      for i := 0; i < p.count; i += 1 {
        p.verts[i] = poly.Points[i]
      }
  }
  return p
}