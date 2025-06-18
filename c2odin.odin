package c2odin

import "core:math"
import "core:math/linalg"

// Define constants. 
GJK_ITERS :: 20
MAX_POLYGON_VERTS :: 8
F32_EPSILON2 : f32 : math.F32_EPSILON * math.F32_EPSILON

// The two main functions. 
Collided :: proc(A, B: ^Shape, ax: ^Matrix = nil, bx: ^Matrix = nil) -> bool {
  switch A.Type {
    case Circle:
      switch(B.Type) {
        case Circle:
          return CircletoCircle(cast(^Circle)A, cast(^Circle)B)
        case AABB:
          return CircletoAABB(cast(^Circle)A, cast(^AABB)B)
        case Capsule:
          return CircletoCapsule(cast(^Circle)A, cast(^Capsule)B)
        case Polygon:
          return CircletoPoly(cast(^Circle)A, cast(^Polygon)B, bx)
        case:
          return false
      }
    case AABB:
      switch(B.Type) {
        case Circle:
          return CircletoAABB(cast(^Circle)B, cast(^AABB)A)
        case AABB:
          return AABBtoAABB(cast(^AABB)A, cast(^AABB)B)
        case Capsule:
          return AABBtoCapsule(cast(^AABB)A, cast(^Capsule)B)
        case Polygon:
          return AABBtoPoly(cast(^AABB)A, cast(^Polygon)B, bx)
        case:
          return false
      }
    case Capsule:
      switch(B.Type) {
        case Circle:
          return CircletoCapsule(cast(^Circle)B, cast(^Capsule)A)
        case AABB:
          return AABBtoCapsule(cast(^AABB)B, cast(^Capsule)A)
        case Capsule:
          return CapsuletoCapsule(cast(^Capsule)A, cast(^Capsule)B)
        case Polygon:
          return CapsuletoPoly(cast(^Capsule)A, cast(^Polygon)B, bx)
        case:
          return false
      }
    case Polygon:
      switch(B.Type) {
        case Circle:
          return CircletoPoly(cast(^Circle)B, cast(^Polygon)A, bx)
        case AABB:
          return AABBtoPoly(cast(^AABB)B, cast(^Polygon)A, bx)
        case Capsule:
          return CapsuletoPoly(cast(^Capsule)B, cast(^Polygon)A, bx)
        case Polygon:
          return PolytoPoly(cast(^Polygon)A, cast(^Polygon)B, ax, bx)
        case:
          return false
      }
  }
  return false
}
Collide :: proc(A, B: ^Shape, ax, bx: ^Matrix, m: ^CollisionManifold) {
  switch A.Type {
    case Circle:
      switch(B.Type) {
        case Circle:
          CircletoCircleManifold(cast(^Circle)A, cast(^Circle)B, m)
        case AABB:
          CircletoAABBManifold(cast(^Circle)A, cast(^AABB)B, m)
        case Capsule:
          CircletoCapsuleManifold(cast(^Circle)A, cast(^Capsule)B, m)
        case Polygon:
          CircletoPolyManifold(cast(^Circle)A, cast(^Polygon)B, bx, m)
      }
    case AABB:
      switch(B.Type) {
        case Circle:
          CircletoAABBManifold(cast(^Circle)B, cast(^AABB)A, m)
        case AABB:
          AABBtoAABBManifold(cast(^AABB)A, cast(^AABB)B, m)
        case Capsule:
          AABBtoCapsuleManifold(cast(^AABB)A, cast(^Capsule)B, m)
        case Polygon:
          AABBtoPolyManifold(cast(^AABB)A, cast(^Polygon)B, bx, m)
      }
    case Capsule:
      switch(B.Type) {
        case Circle:
          CircletoCapsuleManifold(cast(^Circle)B, cast(^Capsule)A, m)
        case AABB:
          AABBtoCapsuleManifold(cast(^AABB)B, cast(^Capsule)A, m)
        case Capsule:
          CapsuletoCapsuleManifold(cast(^Capsule)A, cast(^Capsule)B, m)
        case Polygon:
          CapsuletoPolyManifold(cast(^Capsule)A, cast(^Polygon)B, bx, m)
      }
    case Polygon:
      switch(B.Type) {
        case Circle:
          CircletoPolyManifold(cast(^Circle)B, cast(^Polygon)A, bx, m)
        case AABB:
          AABBtoPolyManifold(cast(^AABB)B, cast(^Polygon)A, bx, m)
        case Capsule:
          CapsuletoPolyManifold(cast(^Capsule)B, cast(^Polygon)A, bx, m)
        case Polygon:
          PolytoPolyManifold(cast(^Polygon)A, cast(^Polygon)B, ax, bx, m)
      }
  }
}
CastRay :: proc(A: ^Ray, B: ^Shape, bx: ^Matrix, out: ^RayCast) -> bool {
  switch(B.Type) {
    case Circle:
      return RaytoCircle(A, cast(^Circle)B, out)
    case AABB:
      return RaytoAABB(A, cast(^AABB)B, out)
    case Capsule:
      return RaytoCapsule(A, cast(^Capsule)B, out)
    case Polygon:
      return RaytoPoly(A, cast(^Polygon)B, bx, out)
  }
  return false // Shouldn't be able to hit this
}
