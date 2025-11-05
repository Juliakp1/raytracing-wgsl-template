fn hit_sphere(center: vec3f, radius: f32, r: ray, record: ptr<function, hit_record>, max: f32)
{
  var oc = r.origin - center;
  var a = dot(r.direction, r.direction);
  var half_b = dot(oc, r.direction);
  var c = dot(oc, oc) - radius * radius;
  var discriminant = half_b * half_b - a * c;

  if (discriminant < 0.0)
  {
    record.hit_anything = false;
    return;
  }

  var sqrtd = sqrt(discriminant);
  var root = (-half_b - sqrtd) / a;
  if (root < RAY_TMIN || root > max)
  {
    root = (-half_b + sqrtd) / a;
    if (root < RAY_TMIN || root > max)
    {
      record.hit_anything = false;
      return;
    }
  }

  record.t = root;
  record.p = ray_at(r, root);
  var newNormal = (record.p - center) / radius;
  record.frontface = dot(r.direction, newNormal) < 0.0;
  record.normal = select(-normalize(newNormal), normalize(newNormal), record.frontface);
  record.hit_anything = true;
}

fn hit_quad(r: ray, Q: vec4f, u: vec4f, v: vec4f, record: ptr<function, hit_record>, max: f32)
{
  var n = cross(u.xyz, v.xyz);
  var normal = normalize(n);
  var D = dot(normal, Q.xyz);
  var w = n / dot(n.xyz, n.xyz);

  var denom = dot(normal, r.direction);
  if (abs(denom) < 0.0001)
  {
    record.hit_anything = false;
    return;
  }

  var t = (D - dot(normal, r.origin)) / denom;
  if (t < RAY_TMIN || t > max)
  {
    record.hit_anything = false;
    return;
  }

  var intersection = ray_at(r, t);
  var planar_hitpt_vector = intersection - Q.xyz;
  var alpha = dot(w, cross(planar_hitpt_vector, v.xyz));
  var beta = dot(w, cross(u.xyz, planar_hitpt_vector));

  if (alpha < 0.0 || alpha > 1.0 || beta < 0.0 || beta > 1.0)
  {
    record.hit_anything = false;
    return;
  }

  if (dot(normal, r.direction) > 0.0)
  {
    record.hit_anything = false;
    return;
  }

  record.t = t;
  record.p = intersection;
  record.normal = normal;
  record.hit_anything = true;
}

fn hit_triangle(r: ray, v0: vec3f, v1: vec3f, v2: vec3f, record: ptr<function, hit_record>, max: f32)
{
  var v1v0 = v1 - v0;
  var v2v0 = v2 - v0;
  var rov0 = r.origin - v0;

  var n = cross(v1v0, v2v0);
  var q = cross(rov0, r.direction);

  var d = 1.0 / dot(r.direction, n);

  var u = d * dot(-q, v2v0);
  var v = d * dot(q, v1v0);
  var t = d * dot(-n, rov0);

  if (u < 0.0 || u > 1.0 || v < 0.0 || (u + v) > 1.0)
  {
    record.hit_anything = false;
    return;
  }

  if (t < RAY_TMIN || t > max)
  {
    record.hit_anything = false;
    return;
  }

  record.t = t;
  record.p = ray_at(r, t);
  record.normal = normalize(n);
  record.hit_anything = true;
}

fn hit_box(r: ray, center: vec3f, rad: vec3f, record: ptr<function, hit_record>, t_max: f32)
{
  var m = 1.0 / r.direction;
  var n = m * (r.origin - center);
  var k = abs(m) * rad;

  var t1 = -n - k;
  var t2 = -n + k;

  var tN = max(max(t1.x, t1.y), t1.z);
  var tF = min(min(t2.x, t2.y), t2.z);

  if (tN > tF || tF < 0.0)
  {
    record.hit_anything = false;
    return;
  }

  var t = tN;
  if (t < RAY_TMIN || t > t_max)
  {
    record.hit_anything = false;
    return;
  }

  record.t = t;
  record.p = ray_at(r, t);
  record.normal = -sign(r.direction) * step(t1.yzx, t1.xyz) * step(t1.zxy, t1.xyz);
  record.hit_anything = true;

  return;
}

fn hit_box_rotated(r: ray, center: vec3f, rad: vec3f, rotation: mat3x3<f32>, record: ptr<function, hit_record>, t_max: f32)
{
  var inv_rotation = transpose(rotation);
  var rotated_origin = inv_rotation * (r.origin - center);
  var rotated_direction = inv_rotation * r.direction;
  var rotated_ray = ray(rotated_origin, rotated_direction);

  hit_box(rotated_ray, vec3f(0.0), rad, record, t_max);

  if (record.hit_anything == true)
  {
    record.normal = rotation * record.normal;
  }
}

fn hit_pyramid(pyramid: pyramid, r: ray, record: ptr<function, hit_record>, max: f32)
{
  var half_size = pyramid.base_size / 2.0;
  var apex = pyramid.base_center.xyz + vec3f(0.0, pyramid.height, 0.0);

  // Base quad
  var base_Q = vec4f(pyramid.base_center.x - half_size, pyramid.base_center.y, pyramid.base_center.z - half_size, 1.0);
  var base_u = vec4f(pyramid.base_size, 0.0, 0.0, 0.0);
  var base_v = vec4f(0.0, 0.0, pyramid.base_size, 0.0);
  hit_quad(r, base_Q, base_u, base_v, record, max);
  if (record.hit_anything == true){return;}

  // Side triangles
  var v0 = pyramid.base_center.xyz + vec3f(-half_size, 0.0, -half_size);
  var v1 = pyramid.base_center.xyz + vec3f( half_size, 0.0, -half_size);
  var v2 = pyramid.base_center.xyz + vec3f( half_size, 0.0,  half_size);
  var v3 = pyramid.base_center.xyz + vec3f(-half_size, 0.0,  half_size);

  hit_triangle(r, v0, v1, apex, record, max);
  if (record.hit_anything == true){return;}

  hit_triangle(r, v1, v2, apex, record, max);
  if (record.hit_anything == true){return;}

  hit_triangle(r, v2, v3, apex, record, max);
  if (record.hit_anything == true){return;}

  hit_triangle(r, v3, v0, apex, record, max);
  if (record.hit_anything == true){return;}
}