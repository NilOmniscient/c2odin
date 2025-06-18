package c2odin

import "core:math/linalg"

Shape :: struct {
  Type: typeid,
  Origin: [2]f32,
}

Circle :: struct {
  using shape: Shape,
  Radius: f32,
}

Capsule :: struct {
  using shape: Shape,
  Radius: f32,
  Origin2: [2]f32,
}

Polygon :: struct{
  using shape: Shape,
  Points: [][2]f32,
  Normals: [][2]f32,
  Count: int,
}

AABB :: struct {
  using shape: Shape,
  Max: [2]f32,
}

CreateAABB :: proc(max: [2]f32) -> AABB {
  return AABB {
    shape = Shape {
      Type = AABB,
      Origin = [2]f32{0, 0},
    },
    Max = max
  }
}
CreateCapsule :: proc(radius: f32, endpoint: [2]f32) -> Capsule {
  return Capsule {
    shape = Shape {
      Type = Capsule,
      Origin = [2]f32{0, 0},
    },
    Radius = radius,
    Origin2 = endpoint,
  }
}
CreateCircle :: proc(radius: f32) -> Circle{
  return Circle {
    shape = Shape {
      Type = Circle,
      Origin = [2]f32{0, 0},
    },
    Radius = radius
  }
}
CreatePolygon :: proc(points: [][2]f32) -> Polygon {
  pcount := len(points)
  hull_count := Hull(points, pcount)
  normals := make([][2]f32, pcount)
  Normals(points[:], normals[:])
  return Polygon {
    shape = Shape {
      Type = Polygon,
      Origin = [2]f32{0, 0}
    },
    Points = points,
    Normals = normals,
    Count = hull_count
  }
}
Hull :: proc(verts: [][2]f32, in_count: int) -> int{
  count := in_count
  if (count <= 2) {
    return 0
  }
  count = min(MAX_POLYGON_VERTS, count)

  right := 0
  xmax : f32 = verts[0].x
  for i := 1; i < count; i += 1 {
    x := verts[i].x
    if (x > xmax) {
      xmax = x
      right = i
    } else if (x == xmax) {
      if (verts[i].y < verts[right].y) {
        right = i
      }
    }
  }

  hull: [MAX_POLYGON_VERTS]int
  out_count := 0
  index := right

  for true {
    hull[out_count] = index
    next := 0
    for i := 1; i < count; i += 1 {
      if (next == index) {
        next = i
        continue;
      }
      e1 := verts[next] - verts[hull[out_count]]
      e2 := verts[i] - verts[hull[out_count]]
      c := linalg.vector_cross2(e1, e2)
      if (c < 0) {
        next = i
      }
      if (c == 0 &&  linalg.vector_dot(e2, e2) > linalg.vector_dot(e1, e1)) {
        next = i
      }
    }

    out_count += 1
    index = next
    if (next == right) {
      break
    }
  }

  hull_verts: [MAX_POLYGON_VERTS][2]f32
  for i := 0; i < out_count; i += 1 {
    hull_verts[i] = verts[hull[i]]
  }
  for i := 0; i < out_count; i += 1 {
    verts[i] = hull_verts[i]
  }
  return out_count
}
Normals :: proc(verts: [][2]f32, norms: [][2]f32) {
  count := len(verts)
  i := 0;
  next := 0
  for i < count {
    next += 1
    a := i
    b := (i + 1) if (i + 1) < count else 0
    e := verts[b] - verts[a]
    to_normalize := [2]f32{e.y, -e.x}
    norms[i] = linalg.vector_normalize(CCW90(e))
    i += 1
  }
}
Inflate :: proc(shape: ^Shape, skin_factor: f32){
  switch shape.Type {
    case Circle:
      c := cast(^Circle)shape
      c.Radius += skin_factor
    case AABB:
      bb := cast(^AABB)shape
      factor := [2]f32{skin_factor, skin_factor}
      bb.Origin -= factor
      bb.Max += factor
    case Capsule:
      c := cast(^Capsule)shape
      c.Radius += skin_factor
    case Polygon:
      poly := cast(^Polygon)shape
      InflatePoly(poly, skin_factor)
  }
}
/* Original Quote:
  Inflating a polytype, idea by Dirk Gregorius ~ 2015. Works in both 2D and 3D.
  Reference: Halfspace intersection with Qhull by Brad Barber
    http://www.geom.uiuc.edu/graphics/pix/Special_Topics/Computational_Geometry/half.html
  
  Algorithm Steps:
  1. Find a points within the input poly
  2. Center this point onto the origin.
  3. Adjust the planes by a skin factor.
  4. Compute the dual vert of each plane. Each plane becomes a vertex.
    c2v dual(c2h plane) {return c2v(plane.n.x / plane.d, plane.n.y / plane.d)}
  5. Compute the convext hull of the dual verts. This is called the dual. 
  6. Compute the dual of the dual, this will be the poly to return.
  7. Translate the poly away from the origin by the center point from step 2.
  8. Return the inflated poly.
*/
InflatePoly :: proc(poly: ^Polygon, skin_factor: f32){
  average := poly.Points[0]
  pcount := poly.Count
  for i := 1; i < pcount; i += 1 {
    average += poly.Points[i]
  }
  average /= cast(f32)pcount

  for i := 0; i < pcount; i += 1 {
    poly.Points[i] -= average
  }

  dual := Dual(poly^, skin_factor)
  poly^ = Dual(dual, 0)
  for i := 0; i < pcount; i += 1 {
    poly.Points[i] += average
  }
}

BBVerts :: proc(bb: ^AABB, out: [][2]f32){
  out[0] = bb.Origin
  out[1] = [2]f32{bb.Max.x, bb.Origin.y}
  out[2] = bb.Max
  out[3] = [2]f32{bb.Origin.x, bb.Max.y}
}

Dual :: proc(poly: Polygon, skin_factor: f32) -> Polygon{
  dual_verts := make([][2]f32, poly.Count)
  for i := 0; i < poly.Count; i += 1 {
    n := poly.Normals[i]
    d := linalg.vector_dot(n, poly.Points[i] + skin_factor)
    if (d == 0) {
      dual_verts[i] = [2]f32{0, 0}
    } else {
      dual_verts[i] = n / d
    }
  }
  dual := CreatePolygon(dual_verts)
  // Create the normals. (No need to recompute the hull, as it's still in CCW order)
  Normals(dual.Points[:], dual.Normals[:])
  return dual
}
