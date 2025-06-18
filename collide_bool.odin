package c2odin

import "core:math"
import "core:math/linalg"

// Collisions with bool output. Faster than manifold, but only good for hit detection. 
CircletoCircle :: proc(a, b: ^Circle) -> bool {
  c := a.Origin - b.Origin
  d2 := linalg.vector_dot(c, c)
  r2 := a.Radius + b.Radius
  r2 *= r2
  return d2 < r2
}
CircletoAABB :: proc(a: ^Circle, b: ^AABB) -> bool {
  l := linalg.clamp(a.Origin, b.Origin, b.Max)
  ab := a.Origin - l
  d2 := linalg.vector_dot(ab, ab)
  r2 := a.Radius*a.Radius
  return d2 < r2
}
AABBtoAABB :: proc(a, b: ^AABB) -> bool {
  d0 := b.Max.x < a.Origin.x
  d1 := a.Max.x < b.Origin.x
  d2 := b.Max.y < a.Origin.y
  d3 := a.Max.y < b.Origin.y
  return !(d0 | d1 | d2 | d3)
}
AABBtoPoint :: proc(a: ^AABB, b: [2]f32) -> bool {
  d0 := b.x < a.Origin.x
  d1 := b.y < a.Origin.x  
  d2 := b.x > a.Max.x
  d3 := b.y > a.Max.y
  return !(d0 | d1 | d2 | d3)
}
CircletoPoint :: proc(a: ^Circle, b: [2]f32) -> bool {
  n := a.Origin - b
  d2 := linalg.vector_dot(n, n)
  r2 := a.Radius * a.Radius
  return d2 < r2
}
CircletoCapsule :: proc(a: ^Circle, b: ^Capsule) -> bool {
  n := b.Origin2 - b.Origin
  ap := a.Origin - b.Origin
  da := linalg.vector_dot(ap, n)
  d2: f32

  if(da < 0) {
    d2 = linalg.vector_dot(ap, ap)
  } else {
    db := linalg.vector_dot(a.Origin - b.Origin2, n)
    if(db < 0) {
      e := ap - (n * (da / linalg.vector_dot(n, n)))
      d2 = linalg.vector_dot(e, e)
    }
  }
  r2 := a.Radius + b.Radius
  r2 *= r2
  return d2 < r2
}

AABBtoCapsule :: proc(a: ^AABB, b: ^Capsule) -> bool{
  return GJK(a, b, nil, nil, nil, nil, true, nil, nil) == 0
}
CapsuletoCapsule :: proc(a, b: ^Capsule) -> bool {
  return GJK(a, b, nil, nil, nil, nil, true, nil, nil) == 0
}
CircletoPoly :: proc(a: ^Circle, b: ^Polygon, bx: ^Matrix) -> bool {
  return GJK(a, b, nil, bx, nil, nil, true, nil, nil) == 0
}
AABBtoPoly :: proc(a: ^AABB, b: ^Polygon, bx: ^Matrix) -> bool {
  return GJK(a, b, nil, bx, nil, nil, true, nil, nil) == 0
}
CapsuletoPoly :: proc(a: ^Capsule, b: ^Polygon, bx: ^Matrix) -> bool {
  return GJK(a, b, nil, bx, nil, nil, true, nil, nil) == 0
}
PolytoPoly :: proc(a, b: ^Polygon, ax, bx: ^Matrix) -> bool {
  return GJK(a, b, ax, bx, nil, nil, true, nil, nil) == 0
}
