package c2odin

import "core:math"
import "core:math/linalg"

// Functions for manipulating vectors, halfspaces, and rays. 

// Rotate the vector 90deg clockwise (Origin top-left)
Skew :: proc(p: [2]f32) -> [2]f32 {
  return [2]f32{
    -p.y, p.x
  }
}
// Rotate the vector 90deg counter clockwise (Origin top-left)
CCW90 :: proc(p: [2]f32) -> [2]f32 {
  return [2]f32{
    p.y, -p.x
  }
}

// Vector functions not necessarily in linalg
vec_min :: proc(a, b: [2]f32) -> [2]f32 {
  return [2]f32{min(a.x, b.x), min(a.y, b.y)}
}
vec_max :: proc(a, b: [2]f32) -> [2]f32 {
  return [2]f32{max(a.x, b.x), max(a.y, b.y)}
}
vec_clamp :: proc(a, lo, hi: [2]f32) -> [2]f32 {
  return vec_max(lo, vec_min(a, hi))
}
vec_abs :: proc(a: [2]f32) -> [2]f32 {
  return [2]f32{abs(a.x), abs(a.y)}
}


/* Helper Functions. 
  Legend:
  Mul = multiply. 
  x = Transformation Matrix
  r = Rotation Vector
  m = Rotation Matrix
  v = Vector
  s = Scalar
  h = HalfSpace
  T = Transpose
*/

// Rotation Procedures
Rot :: proc(rad: f32) -> [2]f32 {
  return [2]f32{
    linalg.cos(rad), linalg.sin(rad)
  }
}
RotIdentity :: proc() -> [2]f32 {
  return [2]f32{1, 0}
}
// This isn't needed, but keeping it to just match the port. 
RotX :: proc(r: [2]f32) -> [2]f32 {
  return r
}
RotY :: proc(r: [2]f32) -> [2]f32 {
  return [2]f32{-r.y, r.x}
}
Mulrv :: proc(r, v: [2]f32) -> [2]f32 {
  return [2]f32{
    r.x * v.x - r.y * v.y, r.y*v.x + r.x*v.y
  }
}
MulrvT :: proc(r, v: [2]f32) -> [2]f32 {
  return [2]f32{
    r.x*v.x + r.y*v.y, r.x*v.y - r.y*v.x
  }
}
// Carbon copies of the Mulrv versions. 
Mulrr :: proc(a, b: [2]f32) -> [2]f32 {
  return Mulrv(a, b)
}
MulrrT :: proc(a, b: [2]f32) -> [2]f32 {
  return MulrvT(a, b)
}

Mulmv :: proc(a: Matrix, b: [2]f32) -> [2]f32 {
  return [2]f32{
    a.x.x * b.x + a.y.x * b.y,
    a.x.y * b.x + a.y.y * b.y
  }
}
MulmvT :: proc(a: Matrix, b: [2]f32) -> [2]f32 {
  return [2]f32{
    a.x.x * b.x + a.x.y * b.y,
    a.y.x * b.x + a.y.y * b.y
  }
}
Mulmm :: proc(a, b: Matrix) -> Matrix {
  return Matrix{
    Mulmv(a, b.x),
    Mulmv(a, b.y),
  }
}
MulmmT :: proc(a, b: Matrix) -> Matrix {
  return Matrix{
    MulmvT(a, b.x),
    MulmvT(a, b.y)
  }
}

// Transform Procedures
TransformIdentity :: proc() -> Matrix {
  return Matrix{
    RotIdentity(),
    [2]f32{0, 0},
  }
}
Mulxv :: proc(a: Matrix, b: [2]f32) -> [2]f32 {
  rv := Mulrv(a.x, b)
  return (rv + a.y)
}
MulxvT :: proc(a: Matrix, b: [2]f32) -> [2]f32 {
  return (MulrvT(a.x, (b - a.y)))
}
Mulxx :: proc(a, b: Matrix) -> Matrix {
  return Matrix{
    Mulrr(a.x, b.x),
    Mulrv(a.x, b.y) + a.y,
  }
}
MulxxT :: proc(a, b: Matrix) -> Matrix {
  return Matrix{
    MulrrT(a.x, b.x),
    MulrvT(a.x, b.y-a.y),
  }
}
Transform :: proc(p: [2]f32, radians: f32) -> Matrix {
  return Matrix{
    Rot(radians),
    p
  }
}

// HalfSpace Procedures
Origin :: proc(h: HalfSpace) -> [2]f32 {
  return h.n * h.d
}
Dist :: proc(h: HalfSpace, p: [2]f32) -> f32 {
  return linalg.vector_dot(h.n, p) - h.d
}
Project :: proc(h: HalfSpace, p: [2]f32) -> [2]f32 {
  return p - (h.n * Dist(h, p))
}
Mulxh :: proc(a: Matrix, b: HalfSpace) -> HalfSpace {
  n := Mulrv(a.x, b.n)
  return HalfSpace{
    n = n,
    d = linalg.vector_dot(
      Mulxv(a, Origin(b)),
      n
    )
  }
}
MulxhT :: proc(a: Matrix, b: HalfSpace) -> HalfSpace {
  n := MulrvT(a.x, b.n)
  return HalfSpace{
    n = n,
    d = linalg.vector_dot(
      MulxvT(a, Origin(b)),
      n
    )
  }
}
Intersect :: proc(a: [2]f32, b: [2]f32, da: f32, db: f32) -> [2]f32 {
  return a+((b-a)*(da / (da - db)))
}
AntiNormalFace :: proc(cap: ^Capsule, p: ^Polygon, x: Matrix, face_out: ^int, n_out: ^[2]f32) {
  sep: f32 = -math.F32_MAX
  index := 0
  n := [2]f32{0, 0}
  pcount := p.Count
  for i := 0; i < pcount; i += 1 {
    h := Mulxh(x, PlaneAt(p, i))
    n0 := -h.n
    s := CapsuleSupport(cap, n0)
    d := HSDist(h, s)
    if (d > sep) {
      sep = d
      index = i
      n = n0
    }
  }
  face_out^ = index
  n_out^ = n
}
Incident :: proc(incident: ^[2][2]f32, ip: ^Polygon, ix: Matrix, rn_in_incident_space: [2]f32) {
  index := 0
  min_dot: f32 = math.F32_MAX

  ipcount := ip.Count
  for i := 0; i < ipcount; i += 1 {
    idot := linalg.vector_dot(rn_in_incident_space, ip.Normals[i])
    if (idot < min_dot) {
      min_dot = idot
      index = i
    }
  }
  incident[0] = Mulxv(ix, ip.Points[index])
  if (ipcount != 1) {
    index += 1
  }
  incident[1] = Mulxv(ix, ip.Points[index])
}
CheckFaces :: proc(A, B: ^Polygon, ax, bx: Matrix, face_index: ^int) -> f32 {
  bina := MulxxT(ax, bx)
  ainb := MulxxT(bx, ax)

  sep: f32 = -math.F32_MAX
  index := 0

  for i := 0; i < A.Count; i += 1 {
    h := PlaneAt(A, i)
    idx := Support(B.Points[:], B.Count, Mulrv(ainb.x, -h.n))
    p := Mulxv(bina, B.Points[idx])
    d := HSDist(h, p)

    if (d > sep) {
      sep = d
      index = i
    }
  }

  face_index^ = index
  return sep
}

// Ray functions. 
Impact :: proc(ray: ^Ray, t: f32) -> [2]f32 {
  return ray.Origin + (ray.Direction * t)
}
SignedDistPointToPlane_1D :: proc(p, n, d: f32) -> f32 {
  return p * n - d * n
}
RayToPlane_1D :: proc(da, db: f32) -> f32 {
  if(da < 0) {
    return 0
  } else if (da * db > 0) {
    return 1.0
  } else {
    d := da - db
    if (d != 0) {
      return da / d
    }
  }
  return 0
}

// Planar functions. 
PlaneAt :: proc(p: ^Polygon, i: int) -> HalfSpace {
  h: HalfSpace
  h.n = p.Normals[i]
  h.d = linalg.vector_dot(p.Normals[i], p.Points[i])
  return h
}
HSDist :: proc(h: HalfSpace, p: [2]f32) -> f32 {
  return linalg.vector_dot(h.n, p) - h.d
}
HSProject :: proc(h: HalfSpace, p: [2]f32) -> [2]f32 {
  return p - (h.n * HSDist(h, p))
}
HSClip :: proc(seg: ^[2][2]f32, h: HalfSpace) -> int {
  out: [2][2]f32
  sp: int = 0
  d0 := HSDist(h, seg[0])
  d1 := HSDist(h, seg[1])
  if (d0 < 0) {
    out[sp] = seg[0]
    sp += 1 
  }
  if (d1 < 0) {
    out[sp] = seg[1]
    sp += 1
  }
  if (d0 == 0 && d1 == 0) {
    out[sp] = seg[0]
    sp += 1
    out[sp] = seg[1]
    sp += 1
  } else if (d0 * d1 <= 0) {
    out[sp] = Intersect(seg[0], seg[1], d0, d1)
    sp += 1
  }
  seg[0] = out[0]
  seg[1] = out[1]
  return sp
}
HSSidePlanes :: proc(seg: ^[2][2]f32, ra, rb: [2]f32, h: ^HalfSpace) -> bool {
  in_vec: [2]f32 = linalg.vector_normalize(rb - ra)
  left := HalfSpace {
    -in_vec,
    linalg.vector_dot(-in_vec, ra)
  }
  right := HalfSpace {
    in_vec,
    linalg.vector_dot(in_vec, rb)
  }
  if (HSClip(seg, left) < 2 || HSClip(seg, right) < 2) {
    return false
  }
  if(h != nil) {
    h.n = CCW90(in_vec)
    h.d = linalg.vector_dot(h.n, ra)
  }
  return true
}
/*
  Per Randy Gaul:
  Clip a segment to the "side planes" of another segment.
  Side Planes are planes orthogonal to a segment, and attached
  to the enpoints fo the segment. 
*/
HSSidePlanesFromPoly :: proc(seg: ^[2][2]f32, x: Matrix, p: ^Polygon, e: int, h: ^HalfSpace) -> bool {
  ie := e
  ra := Mulxv(x, p.Points[ie])
  if(p.Count != 1) {
    ie += 1
  }
  rb := Mulxv(x, p.Points[ie])
  return HSSidePlanes(seg, ra, rb, h)
}
HSKeepDeep :: proc(seg: ^[2][2]f32, h: HalfSpace, m: ^CollisionManifold) {
  cp: u32 = 0
  for i := 0; i < 2; i += 1 {
    p := seg[i]
    d := HSDist(h, p)
    if (d <= 0) {
      m.ContactPoints[cp] = p
      m.Depths[cp] = -d
      cp += 1
    }
    m.Count = cp
    m.Direction = h.n
  }
}
