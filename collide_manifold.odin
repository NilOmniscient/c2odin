package c2odin

import "core:math"
import "core:math/linalg"


CollisionManifold :: struct {
  Count: u32,
  Depths: [2]f32,
  ContactPoints: [2][2]f32,
  Direction: [2]f32,
}
/*
  Collision functions that provide manifolds instead of bools.
  Slower, but much more informative for a proper physics system.
*/
CircletoCircleManifold :: proc(A, B: ^Circle, m: ^CollisionManifold) {
  m.Count = 0
  d := B.Origin - A.Origin
  d2 := linalg.vector_dot(d, d)
  r := A.Radius + B.Radius
  if (d2 < r*r) {
    l := math.sqrt(d2)
    n := d * 1/l if l != 0 else [2]f32{0, 1}
    m.Count = 1
    m.Depths[0] = r - l
    m.ContactPoints[0] = B.Origin - (n*B.Radius)
    m.Direction = n
  }
}

CircletoAABBManifold :: proc(A: ^Circle, B: ^AABB, m: ^CollisionManifold) {
  m.Count = 0
  L: [2]f32 = vec_clamp(A.Origin, B.Origin, B.Max)
  ab := L - A.Origin
  d2 := linalg.vector_dot(ab, ab)
  r2 := A.Radius*A.Radius
  if (d2 < r2) {
    // Shallow case (center of circle not within AABB)
    if (d2 != 0) {
      d := math.sqrt(d2)
      n := linalg.vector_normalize(ab)
      m.Count = 1
      m.Depths[0] = A.Radius - d
      m.ContactPoints[0] = A.Origin + (n*d)
      m.Direction = n
    } else {
      // Deep case (center of circle inside AABB)
      mid := (B.Origin + B.Max) * 0.5
      e := (B.Max - B.Origin) * 0.5
      d := A.Origin - mid
      abs_d := vec_abs(d)
      x_overlap: f32 = e.x - abs_d.x
      y_overlap: f32 = e.y - abs_d.y

      depth: f32
      n: [2]f32

      if (x_overlap < y_overlap) {
        depth = x_overlap
        n = [2]f32{1, 0}
        n *= 1 if d.x < 0 else -1
      } else {
        depth = y_overlap
        n = [2]f32{0, 1}
        n *= 1 if d.y < 0 else -1
      }
      m.Count = 1
      m.Depths[0] = A.Radius + depth
      m.ContactPoints[0] = A.Origin - (n*depth)
      m.Direction = n
    }
  }
}

CircletoCapsuleManifold :: proc(A: ^Circle, B: ^Capsule, m: ^CollisionManifold) {
  m.Count = 0
  a, b: [2]f32
  r: f32 = A.Radius + B.Radius
  d: f32 = GJK(A, B, nil, nil, &a, &b, false, nil, nil)

  if (d < r) {
    n: [2]f32 = linalg.vector_normalize(Skew(B.Origin2 - B.Origin)) if d == 0 else linalg.vector_normalize(b - a)
    m.Count = 1
    m.Depths[0] = r - d
    m.ContactPoints[0] = b - (n * B.Radius)
    m.Direction = n
  }
}

AABBtoAABBManifold :: proc(A, B: ^AABB, m: ^CollisionManifold) {
  m.Count = 0
  mid_a := (A.Origin + A.Max) * 0.5
  mid_b := (B.Origin + B.Max) * 0.5
  eA := vec_abs((A.Max - A.Origin) * 0.5)
  eB := vec_abs((B.Max - B.Origin) * 0.5)
  d := mid_b - mid_a

  dx := eA.x + eB.x - math.abs(d.x)
  dy := eA.y + eB.y - math.abs(d.y)

  if(dx == 0 || dy == 0) {
    return
  }

  depth: f32
  n, p: [2]f32
  if (dx < dy) {
    depth = dx
    p = [2]f32{eA.x, 0}
    if (d.x < 0) {
      n = [2]f32{-1, 0}
      p = mid_a - p
    } else {
      n = [2]f32{1, 0}
      p = mid_a + p
    }
  } else {
    depth = dx
    p = [2]f32{eA.y, 0}
    if(d.y < 0) {
      n = [2]f32{0, -1}
      p = mid_a - p
    } else {
      n = [2]f32{0, 1}
      p = mid_a + p
    }
  }
  m.Count = 1
  m.ContactPoints[0] = p
  m.Depths[0] = depth
  m.Direction = n
}

AABBtoCapsuleManifold :: proc(A: ^AABB, B: ^Capsule, m: ^CollisionManifold) {
  m.Count = 0
  p: Polygon
  BBVerts(A, p.Points[:])
  Normals(p.Points[:], p.Normals[:])
  CapsuletoPolyManifold(B, &p, nil, m)
  m.Direction = -m.Direction
}
CapsuletoCapsuleManifold :: proc(A, B: ^Capsule, m: ^CollisionManifold) {
  m.Count = 0
  a, b: [2]f32
  r := A.Radius + B.Radius
  d := GJK(A, B, nil, nil, &a, &b, false, nil, nil)
  if (d < r) {
    n: [2]f32
    if (d == 0) {
      n = linalg.vector_normalize(Skew(A.Origin2 - A.Origin))
    } else {
      n = linalg.vector_normalize(b - a)
    }
    m.Count = 1
    m.Depths[0] = r - d
    m.ContactPoints[0] = b - (n * B.Radius)
    m.Direction = n
  }
}
CircletoPolyManifold :: proc(A: ^Circle, B: ^Polygon, bx_ptr: ^Matrix, m: ^CollisionManifold) {
  m.Count = 0
  a, b: [2]f32
  d := GJK(A, B, nil, bx_ptr, &a, &b, false, nil, nil)
  if (d != 0) {
    n := b - a
    l := linalg.vector_dot(n, n)
    if (l < A.Radius*A.Radius) {
      l = math.sqrt(l)
      m.Count = 1
      m.ContactPoints[0] = b
      m.Depths[0] = A.Radius - l
      m.Direction = n * (1.0 / l)
    }
  } else {
    bx: Matrix = bx_ptr^ if bx_ptr != nil else TransformIdentity()
    sep: f32 = -math.F32_MAX
    index: int = 0
    local: [2]f32 = MulxvT(bx, A.Origin)

    for i := 0; i < B.Count; i += 1 {
      h := PlaneAt(B, i)
      d = HSDist(h, local)
      if (d > A.Radius) {
        return
      }
      if (d > sep) {
        sep = d
        index = i
      }
    }

    h := PlaneAt(B, index)
    p := HSProject(h, local)
    m.Count = 1
    m.ContactPoints[0] = Mulxv(bx, p)
    m.Depths[0] = A.Radius - sep
    m.Direction = -(Mulrv(bx.r, B.Normals[index]))
  }
}

AABBtoPolyManifold :: proc(A: ^AABB, B: ^Polygon, bx: ^Matrix, m: ^CollisionManifold) {
  m.Count = 0
  p: Polygon
  BBVerts(A, p.Points[:])
  Normals(p.Points[:], p.Normals[:])
  PolytoPolyManifold(&p, B, nil, bx, m)
}

CapsuletoPolyManifold :: proc(A: ^Capsule, B: ^Polygon, bx_ptr: ^Matrix, m: ^CollisionManifold) {
  m.Count = 0
  a, b: [2]f32
  d := GJK(A, B, nil, bx_ptr, &a, &b, false, nil, nil)

  if (d < 1.0e-6) {
    bx := bx_ptr^ if bx_ptr != nil else TransformIdentity()
    AinB: Capsule
    AinB.Origin = MulxvT(bx, A.Origin)
    AinB.Origin2 = MulxvT(bx, A.Origin2)
    ab := linalg.vector_normalize(AinB.Origin - AinB.Origin2)

    // Test the axes of the capsule
    ab_h0, ab_h1: HalfSpace
    ab_h0.n = CCW90(ab)
    ab_h0.d = linalg.vector_dot(AinB.Origin, ab_h0.n)
    v0: int = Support(B.Points[:], B.Count, -ab_h0.n)
    s0 := HSDist(ab_h0, B.Points[v0])

    ab_h1.n = Skew(ab)
    ab_h1.d = linalg.vector_dot(AinB.Origin, ab_h1.n)
    v1 := Support(B.Points[:], B.Count, -ab_h1.n)
    s1 := HSDist(ab_h1, B.Points[v1])

    index := 0
    sep: f32 = -math.F32_MAX
    code := 0
    for i := 0; i < B.Count; i += 1 {
      h := PlaneAt(B, i)
      da := linalg.vector_dot(AinB.Origin, -h.n)
      db := linalg.vector_dot(AinB.Origin2, -h.n)
      d := HSDist(h, AinB.Origin) if da > db else HSDist(h, AinB.Origin2)
      if (d > sep) {
        sep = d
        index = i
      }
    }

    // Track axis of minimum separation. 
    if (s0 > sep) {
      sep = s0
      index = v0
      code = 1
    }
    if (s1 > sep) {
      sep = s1
      index = v1
      code = 2
    }
    switch (code) {
      case 0:
        seg: [2][2]f32 = {A.Origin, A.Origin2}
        h: HalfSpace
        if (!HSSidePlanesFromPoly(&seg, bx, B, index, &h)) {
          return
        }
        HSKeepDeep(&seg, h, m)
        m.Direction = -m.Direction
      case 1:
        incident: [2][2]f32
        Incident(&incident, B, bx, ab_h0.n)
        h: HalfSpace
        if (!HSSidePlanes(&incident, AinB.Origin2, AinB.Origin, &h)) {
          return
        }
        HSKeepDeep(&incident, h, m)
      case 2:
        incident: [2][2]f32
        Incident(&incident, B, bx, ab_h1.n)
        h: HalfSpace
        if (!HSSidePlanes(&incident, AinB.Origin, AinB.Origin2, &h)) {
          return
        }
        HSKeepDeep(&incident, h, m)
      case: // This should never trigger, but just in case...
        return
    }
    for i: u32 = 0; i < m.Count; i += 1 {
      m.Depths[i] += A.Radius
    }
  } else if (d < A.Radius) {
    m.Count = 1
    m.Direction = linalg.vector_normalize(a - b)
    m.ContactPoints[0] = a + (m.Direction*A.Radius)
    m.Depths[0] = A.Radius - d
  }
}
PolytoPolyManifold :: proc(A, B: ^Polygon, ax_ptr, bx_ptr: ^Matrix, m: ^CollisionManifold) {
  m.Count = 0
  ax: Matrix = ax_ptr^ if ax_ptr != nil else TransformIdentity()
  bx: Matrix = bx_ptr^ if bx_ptr != nil else TransformIdentity()
  ea, eb: int

  sa: f32 = CheckFaces(A, B, ax, bx, &ea)
  sb: f32 = CheckFaces(B, A, bx, ax, &eb)

  if (sa >= 0 || sb >= 0) {
    return
  }

  rp := new(Polygon)
  ip := new(Polygon)
  // Make sure to clean these up when finished. 
  defer free(rp)
  defer free(ip)

  rx, ix: Matrix
  re: int
  kRelTol: f32 = 0.95
  kAbsTol: f32 = 0.01

  flip: bool = false
  if (sa * kRelTol > sb + kAbsTol) {
    rp = A
    rx = ax
    ip = B
    ix = bx
    re = ea
    flip = false
  } else {
    rp = B
    rx = bx
    ip = A
    ix = ax 
    re = eb
    flip = true
  }

  incident: [2][2]f32
  Incident(&incident, ip, ix, MulrvT(ix.x, Mulrv(rx.x, rp.Normals[re])))
  rh: HalfSpace
  if (!HSSidePlanesFromPoly(&incident, rx, rp, re, &rh)) {
    return
  }
  HSKeepDeep(&incident, rh, m)
  if (flip) {
    m.Direction = -m.Direction
  }
}