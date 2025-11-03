const THREAD_COUNT = 16;
const RAY_TMIN = 0.0001;
const RAY_TMAX = 100.0;
const PI = 3.1415927f;
const FRAC_1_PI = 0.31830987f;
const FRAC_2_PI = 1.5707964f;

@group(0) @binding(0)  
  var<storage, read_write> fb : array<vec4f>;

@group(0) @binding(1)
  var<storage, read_write> rtfb : array<vec4f>;

@group(1) @binding(0)
  var<storage, read_write> uniforms : array<f32>;

@group(2) @binding(0)
  var<storage, read_write> spheresb : array<sphere>;

@group(2) @binding(1)
  var<storage, read_write> quadsb : array<quad>;

@group(2) @binding(2)
  var<storage, read_write> boxesb : array<box>;

@group(2) @binding(3)
  var<storage, read_write> trianglesb : array<triangle>;

@group(2) @binding(4)
  var<storage, read_write> meshb : array<mesh>;

@group(2) @binding(5)
  var<storage, read_write> pyramidsb : array<pyramid>;

struct ray {
  origin : vec3f,
  direction : vec3f,
};

struct sphere {
  transform : vec4f,
  color : vec4f,
  material : vec4f,
};

struct quad {
  Q : vec4f,
  u : vec4f,
  v : vec4f,
  color : vec4f,
  material : vec4f,
};

struct box {
  center : vec4f,
  radius : vec4f,
  rotation: vec4f,
  color : vec4f,
  material : vec4f,
};

struct triangle {
  v0 : vec4f,
  v1 : vec4f,
  v2 : vec4f,
};

struct mesh {
  transform : vec4f,
  scale : vec4f,
  rotation : vec4f,
  color : vec4f,
  material : vec4f,
  min : vec4f,
  max : vec4f,
  show_bb : f32,
  start : f32,
  end : f32,
};

struct pyramid {
  base_center : vec4f,
  height : f32,
  base_size : f32,
  color : vec4f,
  material : vec4f,
};

struct material_behaviour {
  scatter : bool,
  direction : vec3f,
};

struct camera {
  origin : vec3f,
  lower_left_corner : vec3f,
  horizontal : vec3f,
  vertical : vec3f,
  u : vec3f,
  v : vec3f,
  w : vec3f,
  lens_radius : f32,
};

struct hit_record {
  t : f32,
  p : vec3f,
  normal : vec3f,
  object_color : vec4f,
  object_material : vec4f,
  frontface : bool,
  hit_anything : bool,
};

fn ray_at(r: ray, t: f32) -> vec3f
{
  return r.origin + t * r.direction;
}

fn get_ray(cam: camera, uv: vec2f, rng_state: ptr<function, u32>) -> ray
{
  var rd = cam.lens_radius * rng_next_vec3_in_unit_disk(rng_state);
  var offset = cam.u * rd.x + cam.v * rd.y;
  return ray(cam.origin + offset, normalize(cam.lower_left_corner + uv.x * cam.horizontal + uv.y * cam.vertical - cam.origin - offset));
}

fn get_camera(lookfrom: vec3f, lookat: vec3f, vup: vec3f, vfov: f32, aspect_ratio: f32, aperture: f32, focus_dist: f32) -> camera
{
  var camera = camera();
  camera.lens_radius = aperture / 2.0;

  var theta = degrees_to_radians(vfov);
  var h = tan(theta / 2.0);
  var w = aspect_ratio * h;

  camera.origin = lookfrom;
  camera.w = normalize(lookfrom - lookat);
  camera.u = normalize(cross(vup, camera.w));
  camera.v = cross(camera.u, camera.w);

  camera.lower_left_corner = camera.origin - w * focus_dist * camera.u - h * focus_dist * camera.v - focus_dist * camera.w;
  camera.horizontal = 2.0 * w * focus_dist * camera.u;
  camera.vertical = 2.0 * h * focus_dist * camera.v;

  return camera;
}

fn envoriment_color(direction: vec3f, color1: vec3f, color2: vec3f) -> vec3f
{
  var unit_direction = normalize(direction);
  var t = 0.5 * (unit_direction.y + 1.0);
  var col = (1.0 - t) * color1 + t * color2;

  var sun_direction = normalize(vec3(uniforms[13], uniforms[14], uniforms[15]));
  var sun_color = int_to_rgb(i32(uniforms[17]));
  var sun_intensity = uniforms[16];
  var sun_size = uniforms[18];

  var sun = clamp(dot(sun_direction, unit_direction), 0.0, 1.0);
  col += sun_color * max(0, (pow(sun, sun_size) * sun_intensity));

  return col;
}

// ---------------------------------------------------------------- //

fn rotation_matrix(rot: vec3f) -> mat3x3<f32> {
    let cx = cos(rot.x);
    let sx = sin(rot.x);
    let cy = cos(rot.y);
    let sy = sin(rot.y);
    let cz = cos(rot.z);
    let sz = sin(rot.z);

    let Rx = mat3x3<f32>(
        vec3f(1.0, 0.0, 0.0),
        vec3f(0.0, cx, -sx),
        vec3f(0.0, sx, cx)
    );
    let Ry = mat3x3<f32>(
        vec3f(cy, 0.0, sy),
        vec3f(0.0, 1.0, 0.0),
        vec3f(-sy, 0.0, cy)
    );
    let Rz = mat3x3<f32>(
        vec3f(cz, -sz, 0.0),
        vec3f(sz, cz, 0.0),
        vec3f(0.0, 0.0, 1.0)
    );

    return Rz * Ry * Rx;
}

fn transform_point(pos: vec3f, translation: vec3f, rotation: vec3f, scale: vec3f) -> vec3f {
    let R = rotation_matrix(rotation);
    return (R * (pos * scale)) + translation;
}

// ---------------------------------------------------------------- //

fn check_ray_collision(r: ray, max: f32) -> hit_record
{
  var spheresCount = i32(uniforms[19]);
  var quadsCount = i32(uniforms[20]);
  var boxesCount = i32(uniforms[21]);
  var trianglesCount = i32(uniforms[22]);
  var meshCount = i32(uniforms[27]);
  var pyramidsCount = i32(uniforms[28]);

  var record = hit_record(RAY_TMAX, vec3f(0.0), vec3f(0.0), vec4f(0.0), vec4f(0.0), false, false);
  var closest = record;

  for (var i = 0; i < spheresCount; i = i + 1)
  {
    hit_sphere(spheresb[i].transform.xyz, spheresb[i].transform.w, r, &record, closest.t);
    if (record.hit_anything == true && record.t < closest.t)
    {
      record.object_color = spheresb[i].color;
      record.object_material = spheresb[i].material;
      closest = record;
    }
  }

  for (var i = 0; i < boxesCount; i = i + 1)
  {
    let translation = boxesb[i].center.xyz;
    let rotation = boxesb[i].rotation.xyz;
    let scale = boxesb[i].radius.xyz;
    let transformed_center = transform_point(vec3f(0.0), translation, rotation, vec3f(1.0));
    let transformed_radius = boxesb[i].radius.xyz * scale;    

    hit_box_rotated(r, transformed_center, transformed_radius, rotation_matrix(rotation), &record, closest.t);
    if (record.hit_anything == true && record.t < closest.t)
    {
      record.object_color = boxesb[i].color;
      record.object_material = boxesb[i].material;
      record.frontface = dot(r.direction, record.normal) < 0.0;
      closest = record;
    }
  }

  for (var i = 0; i < quadsCount; i = i + 1)
  {
    hit_quad(r, quadsb[i].Q, quadsb[i].u, quadsb[i].v, &record, closest.t);
    if (record.hit_anything == true && record.t < closest.t)
    {
      record.object_color = quadsb[i].color;
      record.object_material = quadsb[i].material;
      closest = record;
    }
  }

  for (var i = 0; i < meshCount; i = i + 1){

    let translation = meshb[i].transform.xyz;
    let rotation = meshb[i].rotation.xyz;
    let scale = meshb[i].scale.xyz;
    let min = meshb[i].min.xyz;
    let max = meshb[i].max.xyz;

    let transformed_min = transform_point(min, translation, rotation, scale);
    let transformed_max = transform_point(max, translation, rotation, scale);
    hit_box(r, (transformed_min + transformed_max) / 2.0, (transformed_max - transformed_min) / 2.0, &record, closest.t); // bounding box

    if (record.hit_anything == true) {
      for (var j = i32(meshb[i].start); j < i32(meshb[i].end); j = j + 1) {
        var v0 = trianglesb[j].v0.xyz;
        var v1 = trianglesb[j].v1.xyz;
        var v2 = trianglesb[j].v2.xyz;

        v0 = transform_point(v0, meshb[i].transform.xyz, meshb[i].rotation.xyz, meshb[i].scale.xyz);
        v1 = transform_point(v1, meshb[i].transform.xyz, meshb[i].rotation.xyz, meshb[i].scale.xyz);
        v2 = transform_point(v2, meshb[i].transform.xyz, meshb[i].rotation.xyz, meshb[i].scale.xyz);

        hit_triangle(r, v0, v1, v2, &record, closest.t);
        if (record.hit_anything == true && record.t < closest.t)
        {
          record.object_color = meshb[i].color;
          record.object_material = meshb[i].material;
          record.frontface = dot(r.direction, record.normal) < 0.0;
          closest = record;
        }
      }
    }
  }

  // for (var i = 0; i < pyramidsCount; i = i + 1)
  // {
  //   hit_pyramid(pyramidsb[i], r, &record, closest.t);
  //   if (record.hit_anything == true && record.t < closest.t)
  //   {
  //     record.object_color = pyramidsb[i].color;
  //     record.object_material = pyramidsb[i].material;
  //     closest = record;
  //   }
  // }

  return closest;
}

// ---------------------------------------------------------------- //

fn lambertian(normal : vec3f, absorption: f32, random_sphere: vec3f, rng_state: ptr<function, u32>) -> material_behaviour
{
  var scatter_direction = normal + normalize(random_sphere);
  if (length(scatter_direction) < 0.001)
  {
    scatter_direction = normal;
  }
  return material_behaviour(true, scatter_direction);
}

// ---------------------------------------------------------------- //

fn metal(normal : vec3f, direction: vec3f, fuzz: f32, random_sphere: vec3f) -> material_behaviour
{
  var reflected = reflect(direction, normal);
  return material_behaviour(true, reflected + fuzz * random_sphere);
}

// ---------------------------------------------------------------- //

fn schlick(cosine: f32, ref_idx: f32) -> f32 {
  let r0 = ((1.0 - ref_idx) / (1.0 + ref_idx));
  let r0sq = r0 * r0;
  return r0sq + (1.0 - r0sq) * pow(1.0 - cosine, 5.0);
}

fn dielectric(normal : vec3f, r_direction: vec3f, refraction_index: f32, frontface: bool, random_sphere: vec3f, fuzz: f32, rng_state: ptr<function, u32>) -> material_behaviour
{
  let unit_direction = normalize(r_direction);
  let eta = select(refraction_index, 1.0 / refraction_index, frontface);
  let cos_theta = min(dot(-unit_direction, normal), 1.0);
  let sin_theta = sqrt(1.0 - cos_theta * cos_theta);

  // Total internal reflection
  if (eta * sin_theta > 1.0) {
    let reflected = reflect(unit_direction, normal);
    return material_behaviour(true, reflected + fuzz * random_sphere);
  }
    
  // Decide reflection vs refraction using Schlick's approximation
  let reflect_prob = schlick(cos_theta, refraction_index);
  if (rng_next_float(rng_state) < reflect_prob) {
    let reflected = reflect(unit_direction, normal);
    return material_behaviour(true, reflected + fuzz * random_sphere);
  } else {
    let refracted = refract(unit_direction, normal, eta);
    return material_behaviour(true, refracted + fuzz * random_sphere);
  }
}

// ---------------------------------------------------------------- //

fn emmisive(color: vec3f, light: f32) -> material_behaviour
{
    return material_behaviour(false, vec3f(0.0));
}

// ---------------------------------------------------------------- //

fn trace(r: ray, rng_state: ptr<function, u32>) -> vec3f
{
  var maxbounces = i32(uniforms[2]);
  var light = vec3f(0.0);
  var color = vec3f(1.0);
  var r_ = r;
  
  var backgroundcolor1 = int_to_rgb(i32(uniforms[11]));
  var backgroundcolor2 = int_to_rgb(i32(uniforms[12]));
  var behaviour = material_behaviour(true, vec3f(0.0));

  for (var j = 0; j < maxbounces; j = j + 1)
  {
    var record = check_ray_collision(r_, RAY_TMAX);

    if (record.hit_anything == false)
    {
      light = light + color * envoriment_color(r_.direction, backgroundcolor1, backgroundcolor2);
      break;
    }

    var normal = record.normal;
    var frontface = record.frontface;
    var object_material = record.object_material;
    var random_sphere = rng_next_vec3_in_unit_sphere(rng_state);
    
    if (object_material.w > 0.0) // emmisive
    {
      light = light + record.object_color.xyz * color;
      behaviour = emmisive(record.object_color.xyz, object_material.y);
      break;
    }
    else if (object_material.x < 0.0) // dielectric
    {
      behaviour = dielectric(normal, r_.direction, object_material.z, frontface, random_sphere, object_material.y, rng_state);
    }
    else if (object_material.x >= 0.0) // other materials
    {
      var behaviourLambertian = lambertian(normal, object_material.y, random_sphere, rng_state);
      var behaviourMetal = metal(normal, r_.direction, object_material.y, random_sphere);

      var t = object_material.z; // blending factor
      if (rng_next_float(rng_state) < t) {
          behaviour = behaviourMetal;
      } else {
          behaviour = behaviourLambertian;
      }
    }
    
    r_ = ray(record.p, normalize(behaviour.direction));
    r_.origin = r_.origin + r_.direction * RAY_TMIN; // prevent acne
    color = color * record.object_color.xyz;
  }
  
  return saturate(light);
}

// ---------------------------------------------------------------- //

@compute @workgroup_size(THREAD_COUNT, THREAD_COUNT, 1)
fn render(@builtin(global_invocation_id) id : vec3u)
{
    var rez = uniforms[1];
    var time = u32(uniforms[0]);

    // init_rng (random number generator) we pass the pixel position, resolution and frame
    var rng_state = init_rng(vec2(id.x, id.y), vec2(u32(rez)), time);

    // Get uv
    var fragCoord = vec2f(f32(id.x), f32(id.y));
    var uv = (fragCoord + sample_square(&rng_state)) / vec2(rez);

    // Camera
    var lookfrom = vec3(uniforms[7], uniforms[8], uniforms[9]);
    var lookat = vec3(uniforms[23], uniforms[24], uniforms[25]);

    // Get camera
    var cam = get_camera(lookfrom, lookat, vec3(0.0, 1.0, 0.0), uniforms[10], 1.0, uniforms[6], uniforms[5]);
    var samples_per_pixel = i32(uniforms[4]);

    //var color = vec3(rng_next_float(&rng_state), rng_next_float(&rng_state), rng_next_float(&rng_state));
    var color = vec3(0.0);

    // Loop for each sample per pixel
    for (var i = 0; i < samples_per_pixel; i = i + 1)
    {
      var r = get_ray(cam, uv, &rng_state);
      color = color + trace(r, &rng_state);
    }
    color = color / f32(samples_per_pixel);

    var color_out = vec4(linear_to_gamma(color), 1.0);
    var map_fb = mapfb(id.xy, rez);
    
    // Accumulate the color
    var should_accumulate = uniforms[3];
    if (should_accumulate == 0.0){
      rtfb[map_fb] = vec4f(0.0);
    }

    // Set the color to the framebuffer
    rtfb[map_fb] += color_out;
    fb[map_fb] = rtfb[map_fb]/rtfb[map_fb].w;
}