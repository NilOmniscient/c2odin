package c2odin

import "core:math"
import "core:math/linalg"

Ray :: struct {
  using shape: Shape,
  Direction: [2]f32,
  Distance: f32,
}

CreateRay :: proc(direction: [2]f32, distance: f32) -> Ray {
  return Ray {
    shape = Shape {
      Type = Ray,
      Origin = [2]f32{0, 0},
    },
    Direction = direction,
    Distance = distance,
  }
}

RayCast :: struct {
  Time: f32,
  Normal: [2]f32,
}

// Functions for checking ray collisions.
RaytoCircle :: proc(A: ^Ray, B: ^Circle, out: ^RayCast) -> bool {
  p: [2]f32 = B.Origin
  m: [2]f32 = A.Origin - p
  c: f32 = linalg.vector_dot(m, m) - B.Radius*B.Radius
  b: f32 = linalg.vector_dot(m, A.Direction)
  disc: f32 = b * b - c
  if (disc < 0) {
    return false
  }
  t: f32 = -b - math.sqrt(disc)
  if (t >= 0 && t <= A.Distance) {
    out.Time = t
    impact: [2]f32 = Impact(A, t)
    out.Normal = linalg.vector_normalize(impact - p)
    return true
  }
  return false
}
RaytoAABB :: proc(A: ^Ray, B: ^AABB, out: ^RayCast) -> bool {
  p0 := A.Origin
  p1 := Impact(A, A.Distance)
  a_box := AABB{
    Origin = vec_min(p0, p1),
    Max = vec_max(p0, p1),
  }
  // Test B's axes
  if (!AABBtoAABB(&a_box, B)) {
    return false
  }

  // Test the ray's axes along the segment's normal. 
  ab := p1 - p0
  n := Skew(ab)
  abs_n := vec_abs(n)
  half_extents := (B.Max - B.Origin)*0.5
  center_of_b := (B.Max + B.Origin) * 0.5
  d := abs(linalg.vector_dot(n, (p0 - center_of_b))) - linalg.vector_dot(abs_n, half_extents)

  if (d > 0) {
    return false
  }

  // Calculate intermediate values up front. 
  da0 := SignedDistPointToPlane_1D(p0.x, -1, B.Origin.x)
  db0 := SignedDistPointToPlane_1D(p1.x, -1, B.Origin.x)
  da1 := SignedDistPointToPlane_1D(p0.x, 1, B.Max.x)
  db1 := SignedDistPointToPlane_1D(p1.x, 1, B.Max.x)
  da2 := SignedDistPointToPlane_1D(p0.y, -1, B.Origin.y)
  db2 := SignedDistPointToPlane_1D(p1.y, -1, B.Origin.y)
  da3 := SignedDistPointToPlane_1D(p0.y, 1, B.Max.y)
  db3 := SignedDistPointToPlane_1D(p1.y, 1, B.Max.y)

  t0 := RayToPlane_1D(da0, db0)
  t1 := RayToPlane_1D(da1, db1)
  t2 := RayToPlane_1D(da2, db2)
  t3 := RayToPlane_1D(da3, db3)

  // Calculate the hit predicate, no branching.
  hit0 := t0 <= 1
  hit1 := t1 <= 1
  hit2 := t2 <= 1
  hit3 := t3 <= 1

  hit := hit0 | hit1 | hit2 | hit3

  if (hit) {
    t0 *= 1 if hit0 else 0
    t1 *= 1 if hit1 else 0
    t2 *= 1 if hit2 else 0
    t3 *= 1 if hit3 else 0

    if (t0 >= t1 && t0 >= t2 && t0 >= t3) {
      out.Time = t0 * A.Distance
      out.Normal = [2]f32{-1, 0}
    } else if (t1 >= t0 && t1 >= t2 && t1 >= t3) {
      out.Time = t1 * A.Distance
      out.Normal = [2]f32{1, 0}
    } else if (t2 >= t0 && t2 >= t1 && t2 >= t3) {
      out.Time = t2 * A.Distance
      out.Normal = [2]f32{0, -1}
    } else {
      out.Time = t3 * A.Distance
      out.Normal = [2]f32{0, 1}
    }

    return true
  }
  return false
}
RaytoCapsule :: proc(A: ^Ray, B: ^Capsule, out: ^RayCast) -> bool {
  cap_n := B.Origin2 - B.Origin
  norm := linalg.vector_normalize(cap_n)
  M: Matrix = {
    norm,
    CCW90(norm),
  }

  // Rotate both towards the origin along the Y-axis
  yBb := MulmvT(M, cap_n)
  yAp := MulmvT(M, A.Origin - B.Origin)
  yAd := MulmvT(M, A.Direction)
  yAe := yAp + (yAd * A.Distance)

  capsule_bb := AABB{
    Origin = [2]f32{-B.Radius, 0},
    Max = [2]f32{B.Radius, yBb.y},
  }

  out.Normal = norm
  out.Time = 0

  if (AABBtoPoint(&capsule_bb, yAp)) {
    return true
  } else {
    capsule_a := Circle{
      Origin = B.Origin,
      Radius = B.Radius,
    }
    capsule_b := Circle{
      Origin = B.Origin2,
      Radius = B.Radius,
    }
    if (CircletoPoint(&capsule_a, A.Origin) || (CircletoPoint(&capsule_b, A.Origin))) {
      return true
    }
  }

  if (yAe.x * yAp.x < 0 || min(abs(yAe.x), abs(yAp.x)) < B.Radius) {
    Ca := Circle{
      Origin = B.Origin,
      Radius = B.Radius,
    }
    Cb := Circle{
      Origin = B.Origin2,
      Radius = B.Radius,
    }

    if (abs(yAp.x) < B.Radius) {
      if (yAp.y < 0) {
        return RaytoCircle(A, &Ca, out)
      }
      return RaytoCircle(A, &Cb, out)
    } else {
      c := B.Radius if yAp.x > 0 else -B.Radius
      d := yAe.x - yAp.x
      t := (c - yAp.x) / d
      y := yAp.y + ((yAe.y - yAp.y) * t)
      if(y <= 0) {
        return RaytoCircle(A, &Ca, out)
      }
      if(y >= yBb.y) {
        return RaytoCircle(A, &Cb, out)
      } else {
        out.Normal = M.x if c > 0 else Skew(M.y)
        out.Time = t * A.Distance
        return true
      }
    }
  }

  return false
}
RaytoPoly :: proc(A: ^Ray, B: ^Polygon, bx_ptr: ^Matrix, out: ^RayCast) -> bool {
  bx := bx_ptr^ if bx_ptr != nil else TransformIdentity()
  p := MulxvT(bx, A.Origin)
  d := MulrvT(bx.y, A.Direction)
  lo: f32 = 0
  hi: f32 = A.Distance
  index: int = 0

  for i := 0; i < B.Count; i += 1 {
    num: f32 = linalg.vector_dot(B.Normals[i], B.Points[i] - p)
    den: f32 = linalg.vector_dot(B.Normals[i], d)
    if (den == 0 && num < 0) {
      return false
    } else {
      if (den < 0 && num < lo * den) {
        lo = num / den
        index = i
      } else if (den > 0 && num < hi * den) {
        hi = num / den
      }
    }
    if (hi < lo) {
      return false
    }
  }
  if (index != 0) {
    out.Time = lo
    out.Normal = Mulrv(bx.y, B.Normals[index])
    return true
  }
  return false
}
