package c2odin

import "core:math"
import "core:math/linalg"

// Functions relating to GJK and TOI (TOI uses GJK directly)

GJKCache :: struct {
  metric, div: f32,
  count: int,
  iA, iB: [3]int,  
}
GJK :: proc(a, b: ^Shape, ax_ptr, bx_ptr: ^Matrix, out_a, out_b: ^[2]f32, use_radius: bool, iterations: ^int, cache: ^GJKCache) -> f32 {
  ax := ax_ptr^ if ax_ptr != nil else TransformIdentity()
  bx := bx_ptr^ if bx_ptr != nil else TransformIdentity()

  pa := MakeProxy(a)
  pb := MakeProxy(b)

  s := Simplex{}
  verts := &s.verts

  /*
    Explanation via Randy Gaul:
    Metric and caching system as designed by E. Catto in Box2D for his conservative advancement/bilateral
    advancement algorithm implementations. The purpose is to reuse old simplex indices (any simplex that
    has not degenerated into a line or point) as a starting point. This skips the first few iterations of
    GJK going from point, to line, to triangle, thereby lowering the convergence ratees dramatically for
    temporally coherent cases (e.g. time of impact searches).
  */
  cache_was_read := false
  if(cache != nil) {
    // Convert cache count to a bool
    cache_was_good := cache.count > 0
    if(cache_was_good) {
      for i := 0; i < cache.count; i += 0 {
        ia := cache.iA[i]
        ib := cache.iB[i]
        sa := Mulxv(ax, pa.verts[ia])
        sb := Mulxv(bx, pb.verts[ib])
        v: ^SV = &verts[i]
        v.iA = ia
        v.sA = sa
        v.iB = ib
        v.sB = sb
        v.p = sb-sa
        v.u = 0
      }
      s.count = cache.count
      s.div = cache.div

      metric_old := cache.metric
      metric := GJKSimplexMetric(&s)
      min_met := metric if metric < metric_old else metric_old
      max_met := metric if metric > metric_old else metric_old

      if(!(min_met < max_met*2 && metric < -1.0e8)) {
        cache_was_read = true
      }
    }
  }

  if(!cache_was_read) {
    verts[0].iA = 0
    verts[0].iB = 0
    verts[0].sA = Mulxv(ax, pa.verts[0])
    verts[0].sB = Mulxv(bx, pb.verts[0])
    verts[0].p = verts[0].sB - verts[0].sA
    verts[0].u = 1
    s.div = 1
    s.count = 1
  }

  save_a, save_b: [3]int
  save_count: int = 0
  d0: f32 = math.F32_MAX
  d1: f32 = math.F32_MAX
  iter := 0
  hit := false

  for iter < GJK_ITERS {
    save_count = s.count
    for i := 0; i < save_count; i += 1 {
      save_a[i] = verts[i].iA
      save_b[i] = verts[i].iB
    }
    if(s.count == 2){
      c22(&s)
    } else if(s.count == 3){
      c23(&s)
    }
    if(s.count == 3){
      hit = true
      break
    }
    p: [2]f32 = L(&s)
    d1 = linalg.vector_dot(p, p)
    if(d1 > d0) {
      break
    }

    d0 = d1

    d := D(&s)
    if(linalg.vector_dot(d, d) < F32_EPSILON2) {
      break
    }

    ia := Support(pa.verts[:], pa.count, MulrvT(ax.x, -d))
    sa := Mulxv(ax, pa.verts[ia])
    ib := Support(pb.verts[:], pb.count, MulrvT(bx.x, d))
    sb := Mulxv(bx, pb.verts[ib])

    v: ^SV = &verts[s.count]
    v.iA = ia
    v.sA = sa
    v.iB = ib
    v.sB = sb
    v.p = sb-sa

    dup := false
    for i := 0; i < save_count; i += 1 {
      if(ia == save_a[i] && ib == save_b[i]) {
        dup = true
        break
      }
    }
    if(dup) {
      break
    }
    s.count += 1
    iter += 1
  }
  a, b: [2]f32
  Witness(&s, &a, &b)
  dist := linalg.vector_length(a - b)
  if (hit) {
    a = b
    dist = 0
  } else if (use_radius) {
    ra := pa.radius
    rb := pb.radius
    if(dist > (ra + rb) && dist > math.F32_EPSILON) {
      dist -= ra + rb
      n := linalg.vector_normalize(b-a)
      a = a + n*ra
      b = b - n*rb
      if(a.x == b.x && a.y == b.y) {
        dist = 0
      }
    } else {
      p := (a+b)*0.5
      a = p
      b = p
      dist = 0
    }
  }

  if(cache != nil) {
    cache.metric = GJKSimplexMetric(&s)
    cache.count = s.count
    for i := 0; i < s.count; i += 1 {
      v: ^SV = &verts[i]
      cache.iA[i] = v.iA
      cache.iB[i] = v.iB
    }
    cache.div = s.div
  }
  if(out_a != nil) {
    out_a^ = a
  }
  if(out_b != nil) {
    out_b^ = b
  }
  if(iterations != nil) {
    iterations^ = iter
  }
  return dist
}

Witness :: proc(s: ^Simplex, a, b: ^[2]f32) {
  // Original used the Bary macros from the L() func, have to manually expand. 
  den := 1.0 / s.div
  sa := s.verts[s.count].sA
  sb := s.verts[s.count].sB

  switch(s.count) {
    case 1:
      a^ = s.verts[0].sA
      b^ = s.verts[0].sB
    case 2:
      // Bary 2 stuff. 
      a^ = s.verts[0].sA*(den*s.verts[0].u) + s.verts[1].sA*(den*s.verts[1].u)
      b^ = s.verts[0].sB*(den*s.verts[0].u) + s.verts[1].sB*(den*s.verts[1].u)
    case 3:
      // Bary3 stuff. 
      a^ = s.verts[0].sA*(den*s.verts[0].u) + s.verts[1].sA*(den*s.verts[1].u) + s.verts[2].sA*(den*s.verts[2].u)
      b^ = s.verts[0].sB*(den*s.verts[0].u) + s.verts[1].sB*(den*s.verts[1].u) + s.verts[2].sB*(den*s.verts[2].u)
    case:
      a^ = [2]f32{0, 0}
      b^ = [2]f32{0, 0}
  }
  return
}
GJKSimplexMetric :: proc(s: ^Simplex) -> f32 {
  out: f32 = 0
  if(s.count == 2){
    out = linalg.vector_length(s.verts[1].p - s.verts[0].p)
  }
  if(s.count == 3){
    out = linalg.vector_cross2(s.verts[1].p - s.verts[0].p, s.verts[2].p - s.verts[0].p)
  }
  return out
}
c22 :: proc(s: ^Simplex) {
  a: [2]f32 = s.verts[0].p
  b: [2]f32 = s.verts[1].p

  u := linalg.vector_dot(b, b-a)
  v := linalg.vector_dot(a, a-b)

  if v <= 0 {
    s.verts[0].u = 1
    s.div = 1
    s.count = 1
  } else if u <= 0 {
    s.verts[0] = s.verts[1]
    s.verts[0].u = 1
    s.div = 1
    s.count = 1
  } else {
    s.verts[0].u = u
    s.verts[1].u = v
    s.div = u+v
    s.count = 2
  }
  return
}
c23 :: proc(s: ^Simplex) {
  a: [2]f32 = s.verts[0].p
  b: [2]f32 = s.verts[1].p
  c: [2]f32 = s.verts[2].p

  uab := linalg.vector_dot(b, b-a)
  vab := linalg.vector_dot(a, a-b)
  ubc := linalg.vector_dot(c, c-b)
  vbc := linalg.vector_dot(b, b-c)
  uca := linalg.vector_dot(a, a-c)
  vca := linalg.vector_dot(c, c-a)
  area := linalg.vector_cross2(b-a, c-a)
  uabc := linalg.vector_cross2(b, c) * area
  vabc := linalg.vector_cross2(c, a) * area
  wabc := linalg.vector_cross2(a, b) * area

  if (vab <= 0 && uca <= 0) {
    s.verts[0].u = 1
    s.div = 1
    s.count = 1
  } else if (uab <= 0 && vbc <= 0) {
    s.verts[0] = s.verts[1]
    s.verts[0].u = 1
    s.div = 1
    s.count = 1
  } else if (ubc <= 0 && vca <= 0) {
    s.verts[0] = s.verts[2]
    s.verts[0].u = 1
    s.div = 1
    s.count = 1
  } else if (uab > 0 && vab > 0 && wabc <= 0) {
    s.verts[0].u = uab
    s.verts[1].u = vab
    s.div = uab + vab
    s.count = 2
  } else if (ubc > 0 && vbc > 0 && uabc <= 0) {
    s.verts[0] = s.verts[1]
    s.verts[1] = s.verts[2]
    s.verts[0].u = ubc
    s.verts[1].u = vbc
    s.div = ubc + vbc
    s.count = 2
  } else if (uca > 0 && vca > 0 && vabc <= 0) {
    s.verts[1] = s.verts[0]
    s.verts[0] = s.verts[2]
    s.verts[0].u = uca
    s.verts[1].u = vca
    s.div = uca + vca
    s.count = 2
  } else {
    s.verts[0].u = uabc
    s.verts[1].u = vabc
    s.verts[2].u = wabc
    s.div = uabc + vabc + wabc
    s.count = 3
  }
  return
}
L :: proc(s: ^Simplex) -> [2]f32 {
  den := 1.0 / s.div
  if (s.count == 1) {
    return s.verts[0].p
  } else if (s.count == 2) {
    /*
      Original code used a macro:
      BARY(n, x) c2Mulvs(s->n.x, (den*s->n.u)
      BARY2(x) c2Add(BARY(a, x), BARY(b, x))
      BARY3(x) c2Add(c2Add(BARY(a, x), BARY(b, x)), BARY(c, x))
    */
    return (s.verts[0].p * (den * s.verts[0].u)) + (s.verts[1].p * (den * s.verts[1].u))
  }
  return [2]f32{0, 0}
}
D :: proc(s: ^Simplex) -> [2]f32 {
  if(s.count == 1) {
    return -s.verts[0].p
  }
  if(s.count == 2) {
    ab := s.verts[1].p - s.verts[0].p
    if (linalg.vector_cross2(ab, -s.verts[0].p) > 0) {
      return Skew(ab)
    }
    return CCW90(ab)
  }
  return [2]f32{0, 0}
}
TOIResult :: struct {
  hit: bool,
  toi: f32,
  n, p: [2]f32,
  iter: int
}
TOI :: proc(A, B: ^Shape, ax_ptr, bx_ptr: ^Matrix, vA, vB: [2]f32, use_radius: bool) -> TOIResult {
  t: f32 = 0
  ax := ax_ptr^ if ax_ptr != nil else TransformIdentity()
  bx := bx_ptr^ if bx_ptr != nil else TransformIdentity()

  pA := MakeProxy(A)
  pB := MakeProxy(B)

  s: Simplex
  s.count = 0
  verts := &s.verts

  rv := vB - vA
  iA := Support(pA.verts[:], pA.count, MulrvT(ax.x, -rv))
  iB := Support(pB.verts[:], pB.count, MulrvT(bx.x, rv))
  sA := Mulxv(ax, pA.verts[iA])
  sB := Mulxv(bx, pB.verts[iB])
  v := sA - sB

  // Odin has 0 init, so these stay 0 if we aren't using them
  rA, rB, radius: f32
  if(use_radius) {
    rA := pA.radius
    rB := pB.radius
    radius := rA + rB
  }
  tolerance: f32 = 1.0e-4

  result := TOIResult {
    hit = false,
    n = [2]f32{0, 0},
    p = [2]f32{0, 0},
    toi = 1.0,
    iter = 0,
  }

  if(!(linalg.vector_length(v) - radius > tolerance)) {
    result.toi = 0
    result.hit = true
    return result
  }

  for result.iter < 20 {
    iA = Support(pA.verts[:], pA.count, MulrvT(ax.x, -v))
    iB = Support(pB.verts[:], pB.count, MulrvT(bx.x, v))
    sA = Mulxv(ax, pA.verts[iA])
    sB = Mulxv(bx, pB.verts[iB])

    p := sA - sB
    v = linalg.vector_normalize(v)
    vp := linalg.vector_dot(v, p) - radius
    vr := linalg.vector_dot(v, rv)

    if (vp > t * vr) {
      if (vr <= 0) {
        return result
      }
      t = vp / vr
      if (t > 1.0) {
        return result
      }
      result.n = -v
      s.count = 0
    }

    sv: ^SV = &verts[s.count]
    sv.iA = iB
    sv.sA = sB + (rv*t)
    sv.iB = iA
    sv.sB = sA
    sv.p = sv.sB - sv.sA
    sv.u = 1.0
    s.count += 1

    switch(s.count) {
      case 2:
        c22(&s)
      case 3:
        c23(&s)
    }

    if (s.count == 3) {
      result.toi = t
      result.hit = true
      return result
    }
    v = L(&s)
    result.iter += 1
  }

  if (result.iter == 0) {
    result.hit = false
  } else {
    if (linalg.vector_dot(v, v) > 0) {
      result.n = linalg.vector_normalize(-v)
    }
    i := Support(pA.verts[:], pA.count, MulrvT(ax.x, result.n))
    p := Mulxv(ax, pA.verts[i]) + (result.n * rA) + (vA * t)
    result.p = p
    result.toi = t
    result.hit = true
  }
  return result
}