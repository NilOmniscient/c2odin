package c2odin

import "core:math/linalg"

// Support functions
Support :: proc(verts: [][2]f32, count: int, d: [2]f32) -> int {
  imax: int = 0;
  dmax: f32 = linalg.vector_dot(verts[0], d)
  for i := 1; i < count; i += 1 {
    dt := linalg.vector_dot(verts[i], d)
    if(dt > dmax) {
      imax = i
      dmax = dt
    }
  }
  return imax
}
CapsuleSupport :: proc(A: ^Capsule, dir: [2]f32) -> [2]f32 {
  da := linalg.vector_dot(A.Origin, dir)
  db := linalg.vector_dot(A.Origin2, dir)
  if (da > db) {
    return A.Origin + (dir*A.Radius)
  }
  return A.Origin2 + (dir*A.Radius)
}