// Reciprocal, incoherent, opaque-surface geometric rays. NOT full-wave or a film model.
// All lengths are micrometres. MG uses prescribed grey floor and wall/land loss.
// A0 uses the original normalized cross-section curve remapped between endpoints.
// Remapped A0 members are not all claimed to be ideal CPC concentrators.
// Complex electric-field Jones vectors are transported in the global Cartesian frame.
// Convention: exp(-i omega t), passive N=n+i*k, p=s cross propagation_direction.
// Grey active faces use a specified absorptance and PEC-like residual reflection phase.
// Finite trace/cutoff residuals remain an explicit energy bound, never escaped energy.
// [ray_fem_check M3 extension] Surface M3: walls/land are Ti Fresnel; floor ("active") faces
// are a planar coherent single layer (film N1, thickness t) on semi-infinite Ti, evaluated at
// the local incidence cosine per local s/p (Airy). Floor loss is split into film (SiO2) and
// Ti-under-film using the power transmitted into the substrate. Lateral film edges and the
// Ti wall buried inside the film are NOT modelled (planar-stack floor only).
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace {
constexpr double PI=3.1415926535897932384626433832795;
// EPS is conservative bounding-box padding only, not a ray travel distance.
constexpr double EPS=1e-8;
constexpr double INF=std::numeric_limits<double>::infinity();
using Z=std::complex<double>;
struct V {
  double x=0,y=0,z=0;
  V operator+(V b) const {return {x+b.x,y+b.y,z+b.z};}
  V operator-(V b) const {return {x-b.x,y-b.y,z-b.z};}
  V operator*(double a) const {return {a*x,a*y,a*z};}
  double operator[](int i) const {return i==0?x:(i==1?y:z);}
};
double dot(V a,V b){return a.x*b.x+a.y*b.y+a.z*b.z;}
V cross(V a,V b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
double norm(V a){return std::sqrt(dot(a,a));}
V unit(V a){double n=norm(a);if(n==0)throw std::runtime_error("zero vector");return a*(1/n);}
V reflected(V d,V n){return unit(d-n*(2*dot(d,n)));}
struct CV {
  Z x=0,y=0,z=0;
  CV operator+(const CV& b) const {return {x+b.x,y+b.y,z+b.z};}
  CV operator*(double a) const {return {a*x,a*y,a*z};}
};
CV complex_vector(V a,Z b=1){return {a.x*b,a.y*b,a.z*b};}
Z dot(V a,const CV& b){return a.x*b.x+a.y*b.y+a.z*b.z;}
double power(const CV& e){return std::norm(e.x)+std::norm(e.y)+std::norm(e.z);}
double coordinate_scale(V p){return std::max({1.0,std::abs(p.x),std::abs(p.y),std::abs(p.z)});}
double ray_t_min(V p){return 32*std::numeric_limits<double>::epsilon()*coordinate_scale(p);}

struct RNG {
  uint64_t s;
  uint64_t next(){uint64_t z=(s+=UINT64_C(0x9e3779b97f4a7c15));z=(z^(z>>30))*UINT64_C(0xbf58476d1ce4e5b9);z=(z^(z>>27))*UINT64_C(0x94d049bb133111eb);return z^(z>>31);}
  double uniform(){return (next()>>11)*0x1.0p-53;}
};
struct Box {
  V lo={INF,INF,INF},hi={-INF,-INF,-INF};
  void add(V p){lo={std::min(lo.x,p.x),std::min(lo.y,p.y),std::min(lo.z,p.z)};hi={std::max(hi.x,p.x),std::max(hi.y,p.y),std::max(hi.z,p.z)};}
  void add(const Box& b){add(b.lo);add(b.hi);}
  bool hit(V p,V d,double cap) const {
    double a=0,b=cap;
    for(int k=0;k<3;++k){
      if(std::abs(d[k])<1e-15){if(p[k]<lo[k]-EPS||p[k]>hi[k]+EPS)return false;continue;}
      double c=(lo[k]-EPS-p[k])/d[k],e=(hi[k]+EPS-p[k])/d[k];if(c>e)std::swap(c,e);
      a=std::max(a,c);b=std::min(b,e);if(a>b)return false;
    }return b>ray_t_min(p);
  }
};
struct Tri {V a,e1,e2,n,center;Box box;bool active;double area;};
Tri triangle(V a,V b,V c,bool active){
  Tri t;t.a=a;t.e1=b-a;t.e2=c-a;t.n=unit(cross(t.e1,t.e2));t.center=(a+b+c)*(1.0/3);t.active=active;t.area=.5*norm(cross(t.e1,t.e2));t.box.add(a);t.box.add(b);t.box.add(c);return t;
}
double intersect(const Tri& t,V p,V d){
  V q=cross(d,t.e2);double det=dot(t.e1,q);if(std::abs(det)<1e-13)return INF;
  V v=p-t.a;double u=dot(v,q)/det;if(u< -1e-10||u>1+1e-10)return INF;
  V w=cross(v,t.e1);double b=dot(d,w)/det;if(b< -1e-10||u+b>1+1e-10)return INF;
  double dist=dot(t.e2,w)/det;return dist>ray_t_min(p)?dist:INF;
}
double plane_roundoff(const Tri& t,V p){
  double scale=std::max({coordinate_scale(p),coordinate_scale(t.a),coordinate_scale(t.a+t.e1),coordinate_scale(t.a+t.e2)});
  return 64*std::numeric_limits<double>::epsilon()*scale;
}
bool same_reflection_plane(const Tri& a,const Tri& b,V p){
  double tol=std::max(plane_roundoff(a,p),plane_roundoff(b,p));
  return std::abs(dot(a.n,b.n))>=1-64*std::numeric_limits<double>::epsilon()
      &&std::abs(dot(b.a-a.a,a.n))<=tol
      &&std::abs(dot(b.a+b.e1-a.a,a.n))<=tol
      &&std::abs(dot(b.a+b.e2-a.a,a.n))<=tol;
}
bool self_contact(const Tri& candidate,const Tri& previous,V p,V d,double t){
  if(!same_reflection_plane(candidate,previous,p))return false;
  double tol=std::max(plane_roundoff(candidate,p),plane_roundoff(previous,p));
  // A same-plane intersection caused only by the origin's roundoff error:
  // its distance can be amplified at grazing incidence, so use a plane bound.
  return std::abs(dot(p-candidate.a,candidate.n))<=tol
      && t*std::abs(dot(d,candidate.n))<=2*tol;
}
V intersection_point(const Tri& t,V p,V d){
  V q=cross(d,t.e2),v=p-t.a;double det=dot(t.e1,q);
  double u=dot(v,q)/det,b=dot(d,cross(v,t.e1))/det;
  // Reconstruct on the intersected triangle rather than advancing past it.
  return t.a+t.e1*u+t.e2*b;
}
struct Node {Box box;int left=-1,right=-1,first=0,count=0;};
struct Mesh {
  std::vector<Tri> tris;std::vector<int> order;std::vector<Node> nodes;
  void tri(V a,V b,V c,bool active=false){tris.push_back(triangle(a,b,c,active));}
  void quad(V a,V b,V c,V d,bool active=false){tri(a,b,c,active);tri(a,c,d,active);}
  int node(int first,int count){
    int id=static_cast<int>(nodes.size());nodes.emplace_back();Box box,centers;
    for(int i=first;i<first+count;++i){box.add(tris[order[i]].box);centers.add(tris[order[i]].center);}
    nodes[id].box=box;nodes[id].first=first;nodes[id].count=count;
    if(count>8){V ext=centers.hi-centers.lo;int axis=ext.y>ext.x?1:0;if(ext.z>ext[axis])axis=2;
      int mid=first+count/2;std::nth_element(order.begin()+first,order.begin()+mid,order.begin()+first+count,[&](int a,int b){return tris[a].center[axis]<tris[b].center[axis];});
      int l=node(first,count/2),r=node(mid,count-count/2);nodes[id].left=l;nodes[id].right=r;nodes[id].count=0;
    }return id;
  }
  void build(){order.resize(tris.size());std::iota(order.begin(),order.end(),0);nodes.reserve(2*tris.size());if(!tris.empty())node(0,static_cast<int>(tris.size()));}
  int closest(V p,V d,double& distance,int previous_face=-1) const {
    if(nodes.empty())return -1;
    int answer=-1;std::array<int,128> stack{};int n=1;stack[0]=0;
    while(n){const Node& b=nodes[stack[--n]];if(!b.box.hit(p,d,distance))continue;
      if(b.count){for(int j=b.first;j<b.first+b.count;++j){int i=order[j];double t=intersect(tris[i],p,d);if(t<distance){if(previous_face>=0&&self_contact(tris[i],tris[previous_face],p,d,t))continue;distance=t;answer=i;}}}
      else{if(n+2>static_cast<int>(stack.size()))throw std::runtime_error("BVH stack overflow");stack[n++]=b.left;stack[n++]=b.right;}
    }return answer;
  }
  double active_area() const {double a=0;for(const auto&t:tris)if(t.active)a+=t.area;return a;}
};
double cpc_f(double z){double s=std::sin(PI/12),c=std::cos(PI/12),K=20*(1+s);return(-s*(K+c*z)+std::sqrt(K*K+2*K*c*z))/(c*c)-10;}
std::vector<V> polygon(int sides,double radius,double z,double cx=0){
  std::vector<V> p;double phase=sides==6?PI/6:PI/4,circum=radius/std::cos(PI/sides);
  for(int k=0;k<sides;++k){double a=phase+2*PI*k/sides;p.push_back({cx+circum*std::cos(a),circum*std::sin(a),z});}return p;
}
constexpr int ROUND_SIDES=128;
double round_section_area(double vertex_diameter){
  return ROUND_SIDES/8.0*std::sin(2*PI/ROUND_SIDES)*vertex_diameter*vertex_diameter;
}
std::vector<V> round_section(double radius,double z){
  std::vector<V> vertices;
  vertices.reserve(ROUND_SIDES);
  for(int k=0;k<ROUND_SIDES;++k){
    double angle=2*PI*k/ROUND_SIDES;
    vertices.push_back({radius*std::cos(angle),radius*std::sin(angle),z});
  }
  return vertices;
}
std::vector<V> square_section(double radius,double z){
  return {{radius,radius,z},{-radius,radius,z},{-radius,-radius,z},{radius,-radius,z}};
}
struct Geometry {
  std::string name;
  Mesh mesh;
  double height=0,pitch=100,cell_area=10000,opening_area=0,bottom_apothem=10,upper_apothem=35;
  double offset_x=0;
  int opening_sides=0;
  bool cavity=false,hex=false,flat=false,round=false,square=false,trench=false;
  Geometry(const std::string& label,int segments,double requested_height=90,double bottom_width=20,
           double upper_width=2*cpc_f(90),double requested_pitch=80,double requested_offset_x=0):name(label){
    if(!std::isfinite(requested_offset_x))throw std::runtime_error("nonfinite bottom offset");
    if(requested_offset_x!=0&&name!="A0"&&name!="A3")
      throw std::runtime_error("explicit bottom offset is supported only on A0 and A3");
    if(name=="flat_black"||name=="flat_mirror"||name=="flat_patch"||name=="flat_ti"||name=="flat_active") {flat=true;opening_area=cell_area;return;}
    hex=name=="A0"||name=="A1"||name=="A2"||name=="A3";
    round=name=="R0";square=name=="S0";trench=name=="T0";
    if(hex||round||square||trench){
      cavity=true;height=requested_height;pitch=requested_pitch;
      bottom_apothem=bottom_width/2;upper_apothem=upper_width/2;
      if(!std::isfinite(height)||!std::isfinite(pitch)||!std::isfinite(bottom_apothem)||!std::isfinite(upper_apothem)||height<=0||pitch<=0||bottom_apothem<=0||bottom_apothem>=upper_apothem||upper_apothem>=pitch/2)throw std::runtime_error("invalid cavity dimensions or overlapping openings");
      if(std::abs(requested_offset_x)>=upper_apothem-bottom_apothem)
        throw std::runtime_error("bottom offset leaves no positive taper on one side");
      cell_area=(hex||round?std::sqrt(3)*.5:1)*pitch*pitch;
      opening_area=section_area(upper_width);
      opening_sides=hex?6:(round?ROUND_SIDES:(square?4:0));
      offset_x=name=="A1"?-10:(name=="A2"?-20:requested_offset_x);
      if(trench)return; // Open, y-invariant trench is intersected analytically.
      int m=hex&&name!="A3"?segments:1;
      auto previous=section(bottom_apothem,0,offset_x);
      V floor_center={offset_x,0,0};
      for(int k=0;k<opening_sides;++k)
        mesh.tri(floor_center,previous[k],previous[(k+1)%opening_sides],true);
      bool baseline=height==90&&bottom_width==20&&upper_width==2*cpc_f(90);
      double delta=upper_apothem-bottom_apothem;
      for(int j=1;j<=m;++j){
        double z=height*j/m;
        double f=bottom_apothem+delta*z/height;
        if(hex&&name!="A3")
          f=baseline?cpc_f(z):bottom_apothem+delta*(cpc_f(90.0*j/m)-10)/(cpc_f(90)-10);
        double cx=offset_x-offset_x*(f-bottom_apothem)/delta;
        auto current=section(f,z,cx);
        for(int k=0;k<opening_sides;++k){
          int next=(k+1)%opening_sides;
          mesh.quad(previous[k],previous[next],current[next],current[k]);
        }
        previous=current;
      }
    }else if(name=="C"){
      cavity=true;height=40;opening_area=70*70;auto bottom=polygon(4,10,0,-15),top=polygon(4,35,height);V fc={-15,0,0};
      for(int k=0;k<4;++k){mesh.tri(fc,bottom[k],bottom[(k+1)%4],true);mesh.quad(bottom[k],bottom[(k+1)%4],top[(k+1)%4],top[k]);}
    }else if(name=="B"){
      height=40;auto bottom=polygon(4,35,0),top=polygon(4,10,height,-15);V tc={-15,0,height};
      for(int k=0;k<4;++k){mesh.tri(tc,top[k],top[(k+1)%4]);bool plus_x=(bottom[k].x>0&&bottom[(k+1)%4].x>0);mesh.quad(bottom[k],bottom[(k+1)%4],top[(k+1)%4],top[k],plus_x);}
    }else throw std::runtime_error("unknown structure");
    mesh.build();
    if(hex||round||square)validate_mesh();
  }
  std::vector<V> section(double radius,double z,double center_x) const {
    if(round)return round_section(radius,z);
    if(square)return square_section(radius,z);
    return polygon(6,radius,z,center_x);
  }
  double section_area(double width) const {
    if(hex)return std::sqrt(3)*.5*width*width;
    if(round)return round_section_area(width);
    if(trench)return pitch*width;
    return width*width;
  }
  const char* width_definition() const {
    if(hex)return "hex_flat_to_flat";
    if(round)return "polygon_vertex_diameter";
    if(trench)return "trench_clear_x_width";
    if(square)return "square_side";
    return "legacy_geometry";
  }
  const char* profile() const {
    if(name=="A3")return "straight_hexagonal_frustum";
    if(hex)return "cpc_inspired_rescaled_curved_hexagonal_wall";
    if(round)return "straight_regular_128gon_frustum";
    if(square)return "straight_square_frustum";
    if(trench)return "flat_bottom_straight_wall_infinite_y_groove";
    return "legacy_geometry";
  }
  const char* footprint_type() const {
    return hex||round?"triangular_primitive":(trench?"x_periodic_y_invariant_p_by_p_normalization":"square_primitive");
  }
  const char* mesh_validation() const {
    if(hex)return "PASS_manifold_except_six_opening_edges_positive_area_finite_bounds";
    if(round||square)return "PASS_manifold_except_aperture_edges_positive_area_finite_bounds";
    if(trench)return "PASS_analytic_open_y_invariant_xz_sections";
    return "legacy_geometry";
  }
  void validate_mesh() const {
    using Vertex=std::tuple<long long,long long,long long>;
    std::map<std::pair<Vertex,Vertex>,int> edges;
    auto key=[](V p){return Vertex{std::llround(p.x*1e8),std::llround(p.y*1e8),std::llround(p.z*1e8)};};
    for(const auto& t:mesh.tris){
      if(!std::isfinite(t.area)||t.area<=1e-12||std::abs(norm(t.n)-1)>1e-12)throw std::runtime_error("degenerate mesh triangle");
      std::array<V,3> points={t.a,t.a+t.e1,t.a+t.e2};
      for(int j=0;j<3;++j){V p=points[j];if(!std::isfinite(p.x)||!std::isfinite(p.y)||!std::isfinite(p.z)||p.z< -1e-8||p.z>height+1e-8)throw std::runtime_error("mesh vertex outside height range");Vertex a=key(p),b=key(points[(j+1)%3]);if(b<a)std::swap(a,b);++edges[{a,b}];}
    }
    int opening_edges=0;
    for(const auto& item:edges){int count=item.second;if(count==1){const auto& edge=item.first;if(std::abs(std::get<2>(edge.first)*1e-8-height)>1e-7||std::abs(std::get<2>(edge.second)*1e-8-height)>1e-7)throw std::runtime_error("mesh crack away from opening");++opening_edges;}else if(count!=2)throw std::runtime_error("non-manifold mesh edge");}
    if(opening_edges!=opening_sides)throw std::runtime_error("mesh opening topology mismatch");
    double exact=section_area(2*bottom_apothem);
    if(std::abs(active_area()-exact)>1e-9*exact)throw std::runtime_error("mesh floor area mismatch");
  }
  double active_area() const {
    if(trench)return section_area(2*bottom_apothem);
    if(!flat)return mesh.active_area();
    if(name=="flat_black"||name=="flat_active")return cell_area;
    return name=="flat_patch"?2500:0;
  }
  V sample(RNG& r) const {
    double u=r.uniform()-.5,v=r.uniform()-.5;
    return hex||round?V{pitch*u+.5*pitch*v,std::sqrt(3)*.5*pitch*v,height}:V{pitch*u,pitch*v,height};
  }
  bool enter(V& p) const {
    if(trench)return std::abs(p.x)<=upper_apothem;
    if(square)return std::abs(p.x)<=upper_apothem&&std::abs(p.y)<=upper_apothem;
    if(!hex&&!round)return std::abs(p.x)<=35&&std::abs(p.y)<=35;
    double ap=upper_apothem;
    for(int i=-1;i<=1;++i){for(int j=-1;j<=1;++j){V center={pitch*i+.5*pitch*j,std::sqrt(3)*.5*pitch*j,0};V q=p-center;bool inside=true;
      if(round){
        double radius=std::hypot(q.x,q.y),inscribed=ap*std::cos(PI/ROUND_SIDES);
        if(radius>ap+1e-10)inside=false;
        else if(radius>inscribed){
          double sector=2*PI/ROUND_SIDES;
          double delta=std::remainder(std::atan2(q.y,q.x)-sector/2,sector);
          inside=radius*std::cos(delta)<=inscribed+1e-10;
        }
      }else{
        for(int k=0;k<3;++k){
          double a=PI*k/3;
          if(std::abs(q.x*std::cos(a)+q.y*std::sin(a))>ap+1e-10){inside=false;break;}
        }
      }
      if(inside){p=q;return true;}
    }}return false;
  }
};

struct Config {
  std::string structure="A0",surface="M1";
  double theta=0,phi=0,n=0,k=0,grey_alpha=.90,cutoff=1e-10,wavelength=0;
  double height=90,bottom_width=20,upper_width=2*cpc_f(90),pitch=80;
  double offset_x=0,floor_alpha=1,wall_alpha=0;
  // M3 floor film (planar single layer on semi-infinite Ti); lengths in micrometres.
  double film_n=0,film_k=0,film_thickness=0;
  uint64_t rays=10000,seed=1;int segments=128,max_bounces=128;
  bool self_test=false,fresnel_test=false,jones_probe=false,film_probe=false;
  bool floor_alpha_set=false,wall_alpha_set=false;
};
Config args(int argc,char**argv){
  Config c;
  for(int i=1;i<argc;++i){std::string a=argv[i];
    if(a=="--self-test"){c.self_test=true;continue;}if(a=="--fresnel"){c.fresnel_test=true;continue;}if(a=="--jones-probe"){c.jones_probe=true;continue;}if(a=="--film-probe"){c.film_probe=true;continue;}
    if(i+1>=argc)throw std::runtime_error("missing argument value");std::string v=argv[++i];
    if(a=="--structure")c.structure=v;else if(a=="--surface")c.surface=v;else if(a=="--theta")c.theta=std::stod(v);else if(a=="--phi")c.phi=std::stod(v);
    else if(a=="--n")c.n=std::stod(v);else if(a=="--k")c.k=std::stod(v);else if(a=="--grey-alpha")c.grey_alpha=std::stod(v);else if(a=="--cutoff")c.cutoff=std::stod(v);else if(a=="--wavelength-um")c.wavelength=std::stod(v);
    else if(a=="--height-um")c.height=std::stod(v);else if(a=="--bottom-width-um")c.bottom_width=std::stod(v);else if(a=="--upper-width-um")c.upper_width=std::stod(v);else if(a=="--pitch-um")c.pitch=std::stod(v);
    else if(a=="--offset-x-um")c.offset_x=std::stod(v);
    else if(a=="--floor-alpha"){c.floor_alpha=std::stod(v);c.floor_alpha_set=true;}
    else if(a=="--wall-alpha"){c.wall_alpha=std::stod(v);c.wall_alpha_set=true;}
    else if(a=="--film-n")c.film_n=std::stod(v);else if(a=="--film-k")c.film_k=std::stod(v);else if(a=="--film-thickness-um")c.film_thickness=std::stod(v);
    else if(a=="--rays")c.rays=std::stoull(v);else if(a=="--seed")c.seed=std::stoull(v);else if(a=="--segments")c.segments=std::stoi(v);else if(a=="--max-bounces")c.max_bounces=std::stoi(v);else throw std::runtime_error("unknown option: "+a);
  }
  if(c.self_test)return c;
  if(c.film_probe)c.surface="M3";
  if(c.surface!="M0"&&c.surface!="M1"&&c.surface!="M2"&&c.surface!="M3"&&c.surface!="MG")
    throw std::runtime_error("surface must be M0, M1, M2, M3, or MG");
  if(c.surface=="MG"){
    if(!c.floor_alpha_set||!c.wall_alpha_set||!std::isfinite(c.floor_alpha)||!std::isfinite(c.wall_alpha)
       ||c.floor_alpha<0||c.floor_alpha>1||c.wall_alpha<0||c.wall_alpha>1)
      throw std::runtime_error("MG requires explicit --floor-alpha and --wall-alpha in [0,1]");
  }else if(c.floor_alpha_set||c.wall_alpha_set){
    if(c.surface!="M0"||c.floor_alpha!=1||c.wall_alpha!=0)
      throw std::runtime_error("prescribed floor/wall alpha are MG inputs; M0 fixes them at 1/0");
  }
  if(!std::isfinite(c.theta)||!std::isfinite(c.phi)||c.theta<0||c.theta>=90||c.rays==0||c.rays>1000000000ULL||c.segments<1||c.segments>4096||c.max_bounces<1||c.max_bounces>256)throw std::runtime_error("invalid geometry or sampling range");
  if(!std::isfinite(c.grey_alpha)||c.grey_alpha<0||c.grey_alpha>1||!std::isfinite(c.cutoff)||c.cutoff<0||c.cutoff>1e-10||!std::isfinite(c.wavelength)||c.wavelength<0)throw std::runtime_error("invalid alpha, cutoff, or wavelength");
  if(!std::isfinite(c.n)||!std::isfinite(c.k)||c.n<0||c.k<0)throw std::runtime_error("nonfinite or negative refractive index");
  bool need_material=(c.surface!="M0"&&c.surface!="MG")||c.structure=="flat_ti"||c.fresnel_test||c.jones_probe;
  if(need_material&&(!std::isfinite(c.n)||!std::isfinite(c.k)||c.n<=0||c.k<0||(c.k==0&&!c.fresnel_test)))throw std::runtime_error("opaque material requires finite n>0, k>0; --fresnel alone permits k=0");
  if(c.surface=="M3"&&(!std::isfinite(c.film_n)||!std::isfinite(c.film_k)||!std::isfinite(c.film_thickness)||c.film_n<=0||c.film_k<0||c.film_thickness<0||!(c.wavelength>0)))throw std::runtime_error("M3 requires --film-n>0, --film-k>=0, --film-thickness-um>=0 and --wavelength-um>0");
  return c;
}
// Ts,Tp: fraction of the incident local-s/p power transmitted into the semi-infinite substrate
// of a layered face; negative means a single opaque interface (all lost power in that face).
struct Coeff {Z rs,rp;double Ts=-1,Tp=-1;};
Coeff fresnel(double cosine,double n,double k){
  double c=std::max(0.0,std::min(1.0,cosine));Z m(n,k),m2=m*m,q=std::sqrt(m2-(1-c*c));
  // Passive normal wavevector: positive attenuation into the material.
  if(q.imag()<0||(q.imag()==0&&q.real()<0))q=-q;
  return {(c-q)/(c+q),(m2*c-q)/(m2*c+q)};
}
Z passive_q(Z eps,double s2){Z q=std::sqrt(eps-s2);if(q.imag()<0||(q.imag()==0&&q.real()<0))q=-q;return q;}
// Coherent single layer: air (eps0=1) | film eps1=N1^2, thickness d | semi-infinite substrate eps2=N2^2.
// Same sign conventions as fresnel(): s admittance Y=q, r_ij=(Y_i-Y_j)/(Y_i+Y_j);
// p admittance Y=eps/q, r_ij=(Y_j-Y_i)/(Y_j+Y_i) (H-field reflection; reproduces rp of fresnel()).
// r=(r01+r12 e^{2i beta})/(1+r01 r12 e^{2i beta}), beta=(2 pi/lambda) d q1, Im(q_j)>=0.
// Transmitted power into the substrate: s T=Re(q2)/q0 |t_E|^2, t_ij=2q_i/(q_i+q_j);
// p (H form, Z=q/eps) T=Re(q2/eps2)/q0 |t_H|^2, t_ij^H=2 Z_i/(Z_i+Z_j).
Coeff film_stack(double cosine,Z N1,double d_um,Z N2,double lambda_um){
  double c=std::max(0.0,std::min(1.0,cosine));
  if(c<=0){Coeff g{Z(-1),Z(-1)};g.Ts=0;g.Tp=0;return g;} // exactly grazing: total reflection
  double s2=1-c*c;Z e1=N1*N1,e2=N2*N2,q0(c,0),q1=passive_q(e1,s2),q2=passive_q(e2,s2);
  Z beta=(2*PI/lambda_um)*d_um*q1,ph=std::exp(Z(0,1)*beta),ph2=ph*ph;
  Z r01s=(q0-q1)/(q0+q1),r12s=(q1-q2)/(q1+q2),t01s=2.0*q0/(q0+q1),t12s=2.0*q1/(q1+q2);
  Z dens=1.0+r01s*r12s*ph2,rs=(r01s+r12s*ph2)/dens,ts=t01s*t12s*ph/dens;
  Z Y0=1.0/q0,Y1=e1/q1,Y2=e2/q2;
  Z r01p=(Y1-Y0)/(Y1+Y0),r12p=(Y2-Y1)/(Y2+Y1);
  Z Z0=q0,Z1=q1/e1,Z2=q2/e2,t01p=2.0*Z0/(Z0+Z1),t12p=2.0*Z1/(Z1+Z2);
  Z denp=1.0+r01p*r12p*ph2,rp=(r01p+r12p*ph2)/denp,tp=t01p*t12p*ph/denp;
  Coeff f{rs,rp};f.Ts=q2.real()/c*std::norm(ts);f.Tp=Z2.real()/c*std::norm(tp);
  return f;
}
Coeff film_coefficients(const Config& c,double cosine){return film_stack(cosine,Z(c.film_n,c.film_k),c.film_thickness,Z(c.n,c.k),c.wavelength);}
Coeff coefficients(const Config& c,bool active,double cosine,int force=0){
  // force 1=black, 2=perfect mirror, 3=Ti, used only by flat controls.
  if(force==1||(force==0&&c.surface=="M0"&&active))return {0,0};
  if(force==2||(force==0&&c.surface=="M0"))return {-1,1};
  if(force==0&&c.surface=="MG"){
    double amplitude=std::sqrt(1-(active?c.floor_alpha:c.wall_alpha));
    return {-amplitude,amplitude};
  }
  if(force==0&&c.surface=="M2"&&active){double r=std::sqrt(1-c.grey_alpha);return {-r,r};}
  if(force==0&&c.surface=="M3"&&active)return film_coefficients(c,cosine);
  return fresnel(cosine,c.n,c.k);
}
struct Transform {
  V s,p_in,p_out,d_out;Coeff f;
  Transform(V d,V normal,const Config& c,bool active,int force=0){
    V normal_in=dot(d,normal)>0?normal*(-1):normal;
    d_out=reflected(d,normal_in);s=cross(d,normal_in);
    if(norm(s)<1e-12){V axis=std::abs(d.x)<=std::abs(d.y)&&std::abs(d.x)<=std::abs(d.z)?V{1,0,0}:(std::abs(d.y)<=std::abs(d.z)?V{0,1,0}:V{0,0,1});s=cross(d,axis);}
    s=unit(s);p_in=unit(cross(s,d));p_out=unit(cross(s,d_out));f=coefficients(c,active,-dot(d,normal_in),force);
  }
  CV apply(const CV& e) const {return complex_vector(s,f.rs*dot(s,e))+complex_vector(p_out,f.rp*dot(p_in,e));}
};
// FILM and TIFLOOR are sub-tallies of ACTIVE (layered M3 floor: film volume loss and
// Ti-under-film loss); LAND is the sub-tally of WALL deposited on the flat top land.
enum Metric {TOTAL=0,ACTIVE=1,WALL=2,ESCAPE=3,RESIDUAL=4,FILM=5,TIFLOOR=6,LAND=7,SIDEWALL=8,METRICS=9};
using Energy=std::array<double,METRICS>;
double max_split_error=0; // |(film+Ti_floor)-(before-remaining)| over all layered events
struct Result {
  std::array<Energy,2> energy{};std::array<CV,2> e;
  int bounces=0,reason=0;uint64_t wraps=0;
  // reason: 0 exhausted/escaped, 1 intensity cutoff, 2 bounce/event cap, 3 geometry unresolved.
  void finish(bool escape,int why){reason=why;for(int j=0;j<2;++j){energy[j][escape?ESCAPE:RESIDUAL]+=power(e[j]);energy[j][TOTAL]=energy[j][ACTIVE]+energy[j][WALL];}}
};
void interact(Result& r,V& d,V normal,bool active,const Config& c,int force=0,bool land=false){
  Transform tr(d,normal,c,active,force);
  for(int j=0;j<2;++j){double before=power(r.e[j]);CV after=tr.apply(r.e[j]);double remaining=power(after);
    if(!std::isfinite(remaining)||remaining>before+1e-12)throw std::runtime_error("non-passive or nonfinite reflection");
    if(remaining>before){after=after*std::sqrt(before/remaining);remaining=before;}
    double lost=before-remaining;
    r.energy[j][active?ACTIVE:WALL]+=lost;
    if(!active)r.energy[j][land?LAND:SIDEWALL]+=lost;
    if(active&&tr.f.Ts>=0){
      // Planar stack: local s and p are orthogonal eigen-polarizations, so the power entering
      // the substrate is |E_s|^2 Ts + |E_p|^2 Tp; the film takes the remainder of the loss.
      double es=std::norm(dot(tr.s,r.e[j])),ep=std::norm(dot(tr.p_in,r.e[j]));
      double ti=es*tr.f.Ts+ep*tr.f.Tp;
      double film=es*(1-std::norm(tr.f.rs)-tr.f.Ts)+ep*(1-std::norm(tr.f.rp)-tr.f.Tp);
      max_split_error=std::max(max_split_error,std::abs(film+ti-lost));
      if(film< -1e-12||ti< -1e-12)throw std::runtime_error("negative layered absorption");
      r.energy[j][TIFLOOR]+=ti;r.energy[j][FILM]+=lost-ti; // exact bookkeeping: FILM+TIFLOOR==ACTIVE
    }
    r.e[j]=after;
  }
  d=tr.d_out;++r.bounces;
}
void advance_trench(Result& result,V& position,V direction,double distance,double pitch){
  position=position+direction*distance;
  // y is periodic bookkeeping only. It never limits an intersection.
  double periods=std::floor(position.y/pitch+.5);
  position.y-=periods*pitch;
  result.wraps+=static_cast<uint64_t>(std::abs(periods));
}
Result trace_trench(const Geometry& g,V p,V d,const std::array<CV,2>& initial,const Config& c){
  Result r;
  r.e=initial;
  if(!g.enter(p)){
    interact(r,d,{0,0,1},false,c,0,true);
    r.finish(true,0);
    return r;
  }
  double slope=(g.upper_apothem-g.bottom_apothem)/g.height;
  for(;;){
    double ws=power(r.e[0]),wp=power(r.e[1]);
    if(ws==0&&wp==0){r.finish(false,0);return r;}
    if(ws<=c.cutoff&&wp<=c.cutoff){r.finish(false,1);return r;}
    if(r.bounces>=c.max_bounces){r.finish(false,2);return r;}

    double nearest=INF;
    enum class Event {Unresolved,Escape,Floor,Wall};
    Event event=Event::Unresolved;
    int wall_sign=0;
    double tolerance=ray_t_min(p);
    if(d.z>1e-15){
      nearest=std::max(0.0,(g.height-p.z)/d.z);
      event=Event::Escape;
    }else if(d.z< -1e-15){
      double distance=-p.z/d.z;
      double x=p.x+distance*d.x;
      if(distance>=-tolerance&&std::abs(x)<=g.bottom_apothem+tolerance){
        nearest=std::max(0.0,distance);
        event=Event::Floor;
      }
    }
    for(int sign:{-1,1}){
      double outward_speed=sign*d.x-slope*d.z;
      if(outward_speed<=1e-15)continue;
      double distance=(g.bottom_apothem+slope*p.z-sign*p.x)/outward_speed;
      if(distance< -tolerance||distance>=nearest)continue;
      double z=p.z+distance*d.z;
      if(z< -tolerance||z>g.height+tolerance)continue;
      nearest=std::max(0.0,distance);
      event=Event::Wall;
      wall_sign=sign;
    }
    if(event==Event::Unresolved){r.finish(false,3);return r;}
    advance_trench(r,p,d,nearest,g.pitch);
    if(event==Event::Escape){r.finish(true,0);return r;}
    if(event==Event::Floor){
      p.z=0;
      interact(r,d,{0,0,1},true,c);
    }else{
      p.z=std::max(0.0,std::min(g.height,p.z));
      p.x=wall_sign*(g.bottom_apothem+slope*p.z);
      interact(r,d,unit({-double(wall_sign),0,slope}),false,c);
    }
  }
}
Result trace(const Geometry& g,V p,V d,const std::array<CV,2>& initial,const Config& c){
  if(g.trench)return trace_trench(g,p,d,initial,c);
  Result r;r.e=initial;
  if(g.flat){bool active=g.name=="flat_black"||g.name=="flat_active"||(g.name=="flat_patch"&&std::abs(p.x)<=25&&std::abs(p.y)<=25);
    int force=g.name=="flat_black"?1:(g.name=="flat_mirror"?2:(g.name=="flat_ti"?3:0));
    interact(r,d,{0,0,1},active,c,force);r.finish(true,0);return r;
  }
  if(g.cavity){if(!g.enter(p)){interact(r,d,{0,0,1},false,c,0,true);r.finish(true,0);return r;}}
  else p.z+=4*ray_t_min(p); // Begin strictly above protrusion tops, outside solid.
  int previous_face=-1;
  uint64_t event_cap=static_cast<uint64_t>(c.max_bounces)*100+1000;
  for(uint64_t event=0;event<event_cap;++event){
    double ws=power(r.e[0]),wp=power(r.e[1]);
    if(ws==0&&wp==0){r.finish(false,0);return r;}
    if(ws<=c.cutoff&&wp<=c.cutoff){r.finish(false,1);return r;}
    if(r.bounces>=c.max_bounces){r.finish(false,2);return r;}
    double t_top=INF;if(d.z>1e-15)t_top=(g.height-p.z)/d.z;
    if(t_top<0&&d.z>0){r.finish(true,0);return r;}
    double nearest=t_top;int face=g.mesh.closest(p,d,nearest,previous_face);
    double t_floor=INF,t_x=INF,t_y=INF;
    if(!g.cavity){
      if(d.z< -1e-15)t_floor=-p.z/d.z;
      if(d.x>1e-15)t_x=(g.pitch*.5-p.x)/d.x;else if(d.x< -1e-15)t_x=(-g.pitch*.5-p.x)/d.x;
      if(d.y>1e-15)t_y=(g.pitch*.5-p.y)/d.y;else if(d.y< -1e-15)t_y=(-g.pitch*.5-p.y)/d.y;
      double seam=std::min(t_x,t_y);
      if(seam<nearest&&seam<t_floor){
        double tolerance=std::max(ray_t_min(p),32*std::numeric_limits<double>::epsilon()*std::max(1.0,std::abs(seam)));
        if(seam< -tolerance){r.finish(false,3);return r;}p=p+d*std::max(0.0,seam);
        if(std::abs(t_x-seam)<=tolerance)p.x=-std::copysign(g.pitch*.5,d.x);
        if(std::abs(t_y-seam)<=tolerance)p.y=-std::copysign(g.pitch*.5,d.y);
        previous_face=-1; // A translated neighbouring cell has distinct faces.
        ++r.wraps;continue;
      }
      if(t_floor<nearest){if(t_floor< -ray_t_min(p)){r.finish(false,3);return r;}p=p+d*std::max(0.0,t_floor);p.z=0;interact(r,d,{0,0,1},false,c);previous_face=-1;continue;}
    }
    if(face>=0){const Tri&t=g.mesh.tris[face];p=intersection_point(t,p,d);interact(r,d,t.n,t.active,c);previous_face=face;continue;}
    if(std::isfinite(t_top)&&t_top>=0){r.finish(true,0);return r;}
    r.finish(false,3);return r;
  }
  r.finish(false,2);return r;
}
struct Stats {
  uint64_t count=0;std::array<long double,METRICS> mean{},m2{};
  void add(const Energy& x){++count;for(int k=0;k<METRICS;++k){long double delta=x[k]-mean[k];mean[k]+=delta/count;m2[k]+=delta*(x[k]-mean[k]);}}
  double se(int k) const {return count>1?std::sqrt(static_cast<double>(std::max(0.0L,m2[k])/(count-1)/count)):0;}
};
const char* names[METRICS]={"A_total","A_active","A_wall","escape","residual","A_film","A_Ti_floor","A_land","A_sidewall"};
void print_vec(V a){std::cout<<'['<<a.x<<','<<a.y<<','<<a.z<<']';}
void print_complex(Z a){std::cout<<'['<<a.real()<<','<<a.imag()<<']';}
void print_vec(const CV& a){std::cout<<'[';print_complex(a.x);std::cout<<',';print_complex(a.y);std::cout<<',';print_complex(a.z);std::cout<<']';}
void print_stats(const Stats& s){
  std::cout<<'{';for(int k=0;k<METRICS;++k){if(k)std::cout<<',';std::cout<<'"'<<names[k]<<"\":"<<s.mean[k];}
  std::cout<<",\"absorption_upper_bound\":"<<s.mean[TOTAL]+s.mean[RESIDUAL]<<",\"closure\":"<<s.mean[TOTAL]+s.mean[ESCAPE]+s.mean[RESIDUAL]<<",\"se\":{";
  for(int k=0;k<METRICS;++k){if(k)std::cout<<',';std::cout<<'"'<<names[k]<<"\":"<<s.se(k);}std::cout<<"}}";
}
void fresnel_output(const Config& c){Coeff f=fresnel(std::cos(c.theta*PI/180),c.n,c.k);std::cout<<"{\"theta_deg\":"<<c.theta<<",\"n\":"<<c.n<<",\"k\":"<<c.k<<",\"rs\":";print_complex(f.rs);std::cout<<",\"rp\":";print_complex(f.rp);std::cout<<",\"Rs\":"<<std::norm(f.rs)<<",\"Rp\":"<<std::norm(f.rp)<<",\"p_basis\":\"s_cross_d\"}\n";}
void film_output(const Config& c){
  double cs=std::cos(c.theta*PI/180);Coeff f=film_coefficients(c,cs),b=fresnel(cs,c.n,c.k);
  std::cout<<"{\"theta_deg\":"<<c.theta<<",\"wavelength_um\":"<<c.wavelength<<",\"film_n\":"<<c.film_n<<",\"film_k\":"<<c.film_k<<",\"film_thickness_um\":"<<c.film_thickness<<",\"n\":"<<c.n<<",\"k\":"<<c.k<<",\"rs\":";print_complex(f.rs);std::cout<<",\"rp\":";print_complex(f.rp);
  std::cout<<",\"Rs\":"<<std::norm(f.rs)<<",\"Rp\":"<<std::norm(f.rp)<<",\"Ts_substrate\":"<<f.Ts<<",\"Tp_substrate\":"<<f.Tp<<",\"As_film\":"<<1-std::norm(f.rs)-f.Ts<<",\"Ap_film\":"<<1-std::norm(f.rp)-f.Tp;
  std::cout<<",\"bare_rs\":";print_complex(b.rs);std::cout<<",\"bare_rp\":";print_complex(b.rp);std::cout<<",\"p_basis\":\"s_cross_d\",\"p_reflection\":\"H_field_equal_to_fresnel_rp_convention\"}\n";
}
void jones_output(const Config& c){
  V d=unit({.3,-.4,-.8});V s=unit(cross(d,{0,0,1})),p=unit(cross(s,d));
  CV e=complex_vector(s,std::sqrt(.6))+complex_vector(p,Z(.3,.4)*std::sqrt(1.6));
  std::array<V,2> normals={unit({.1,.2,1}),unit({.7,-.3,-1})};
  std::cout<<"{\"n\":"<<c.n<<",\"k\":"<<c.k<<",\"initial_direction\":";print_vec(d);std::cout<<",\"initial_E\":";print_vec(e);std::cout<<",\"normals\":[";print_vec(normals[0]);std::cout<<',';print_vec(normals[1]);std::cout<<"],\"steps\":[";
  for(int j=0;j<2;++j){if(j)std::cout<<',';double before=power(e);Transform tr(d,normals[j],c,false,3);e=tr.apply(e);d=tr.d_out;std::cout<<"{\"direction\":";print_vec(d);std::cout<<",\"E\":";print_vec(e);std::cout<<",\"power_before\":"<<before<<",\"power_after\":"<<power(e)<<'}';}std::cout<<"]}\n";
}
std::array<CV,2> initial_fields(double phi,V d){V s={-std::sin(phi),std::cos(phi),0};return {complex_vector(s),complex_vector(unit(cross(s,d)))};}
void test_new_geometries(){
  Geometry curved("A0",128,82,25,66,80);
  Geometry straight("A3",1,82,25,66,80);
  Geometry round("R0",1,82,25,66,80);
  Geometry square("S0",1,82,25,66,80);
  Geometry trench("T0",1,82,25,66,80);
  Geometry positive("A0",128,82,25,66,80,5);
  Geometry negative("A0",128,82,25,66,80,-5);
  Geometry scaled("A0",128,123,37.5,99,120,7.5);
  Config config;
  config.surface="MG";config.floor_alpha=.9;config.wall_alpha=.04;
  V normal{0,0,-1};
  auto normal_fields=initial_fields(0,normal);
  RNG rng{34576};
  auto verify_energy=[](const Result& result){
    if(result.reason==3)throw std::runtime_error("geometry unresolved in family self-test");
    for(const auto& energy:result.energy){
      if(std::abs(energy[TOTAL]+energy[ESCAPE]+energy[RESIDUAL]-1)>1e-11
         ||std::abs(energy[WALL]-energy[LAND]-energy[SIDEWALL])>1e-12
         ||std::abs(energy[TOTAL]-energy[ACTIVE]-energy[SIDEWALL]-energy[LAND])>1e-12)
        throw std::runtime_error("family energy closure or partition failed");
    }
  };
  for(const Geometry* geometry:{&curved,&straight,&round,&square,&trench}){
    auto direct=trace(*geometry,{0,0,geometry->height},normal,normal_fields,config);
    verify_energy(direct);
    for(const auto& energy:direct.energy)
      if(std::abs(energy[ACTIVE]-.9)>1e-13||std::abs(energy[ESCAPE]-.1)>1e-13)
        throw std::runtime_error("direct normal floor catch failed");
  }
  V direction=unit({-.42,.23,-.8}),mirror_direction={-direction.x,direction.y,direction.z};
  auto fields=initial_fields(std::atan2(-direction.y,-direction.x),direction);
  auto mirror_fields=initial_fields(std::atan2(-mirror_direction.y,-mirror_direction.x),mirror_direction);
  for(const char* surface:{"M0","MG"}){
    config.surface=surface;
    for(const Geometry* geometry:{&curved,&straight,&round,&square,&trench}){
      for(int i=0;i<128;++i)
        verify_energy(trace(*geometry,geometry->sample(rng),direction,fields,config));
    }
    for(int i=0;i<256;++i){
      V p=positive.sample(rng);
      auto a=trace(positive,p,direction,fields,config);
      auto b=trace(negative,{-p.x,p.y,p.z},mirror_direction,mirror_fields,config);
      auto c=trace(scaled,p*1.5,direction,fields,config);
      verify_energy(a);verify_energy(b);verify_energy(c);
      // Mirrored rays and uniformly scaled geometry have the same ray-path absorption.
      for(int polarization=0;polarization<2;++polarization){
        for(int metric:{TOTAL,ACTIVE,WALL,ESCAPE,RESIDUAL}){
          if(std::abs(a.energy[polarization][metric]-b.energy[polarization][metric])>1e-11
             ||std::abs(a.energy[polarization][metric]-c.energy[polarization][metric])>1e-11)
            throw std::runtime_error("mirrored offset or uniform scaling failed");
        }
      }
    }
  }
  auto trench_ray=trace(trench,{0,0,trench.height},unit({0,.99,-.1}),initial_fields(-PI/2,unit({0,.99,-.1})),config);
  verify_energy(trench_ray);
  if(trench_ray.wraps<2||std::abs(trench_ray.energy[0][ACTIVE]-.9)>1e-13)
    throw std::runtime_error("infinite groove y transport failed");
}
void self_test(){
  Config c;c.n=2;c.k=3;c.surface="M1";V d={0,0,-1};auto initial=initial_fields(0,d);Geometry flat("flat_ti",1);
  auto out=trace(flat,{0,0,0},d,initial,c);double exact=1-std::norm((Z(1)-Z(2,3))/(Z(1)+Z(2,3)));
  if(std::abs(out.energy[0][TOTAL]-exact)>1e-13||std::abs(out.energy[1][TOTAL]-exact)>1e-13)throw std::runtime_error("normal Fresnel failed");
  Transform normal(d,{0,0,1},c,false);Coeff f=fresnel(1,2,3);CV ex=normal.apply(complex_vector({1,0,0}));if(std::abs(ex.x-f.rs)>1e-13)throw std::runtime_error("normal global polarization failed");
  Geometry a("A0",128),b("B",1),cavity("C",1);
  if(std::abs(a.active_area()-200*std::sqrt(3))>1e-8||std::abs(b.active_area()-1800*std::sqrt(2))>1e-8||std::abs(cavity.active_area()-400)>1e-8)throw std::runtime_error("geometry areas failed");
  RNG rng{7};for(const char* surface:{"M0","M1","M2"}){c.surface=surface;for(const Geometry* g:{&a,&b,&cavity}){for(int j=0;j<100;++j){auto r=trace(*g,g->sample(rng),d,initial,c);for(const auto& en:r.energy)if(std::abs(en[TOTAL]+en[ESCAPE]+en[RESIDUAL]-1)>1e-11)throw std::runtime_error("energy closure failed");}}}
  c.surface="M2";Geometry active("flat_active",1);auto grey=trace(active,{0,0,0},d,initial,c);if(std::abs(grey.energy[0][ACTIVE]-.9)>1e-13)throw std::runtime_error("grey active failed");
  // --- M3 layered-floor checks (ray_fem_check extension) ---
  // (a) zero thickness reproduces the bare Fresnel coefficients and puts all loss in the substrate.
  for(double th:{0.0,15.0,30.0,45.0,60.0,75.0,85.0,89.0}){double cs=std::cos(th*PI/180);
    for(Z N1:{Z(2.2447212969283283,0.06982232081911283),Z(0.5947997319778178,1.7674602002464557),Z(1,0)}){
      Coeff f=film_stack(cs,N1,0.0,Z(4.308219763636363,20.63042690909091),10.5),b=fresnel(cs,4.308219763636363,20.63042690909091);
      if(std::abs(f.rs-b.rs)>1e-12||std::abs(f.rp-b.rp)>1e-12||std::abs(f.Ts-(1-std::norm(b.rs)))>1e-12||std::abs(f.Tp-(1-std::norm(b.rp)))>1e-12)throw std::runtime_error("zero-thickness film does not reproduce bare Fresnel");}}
  // (b) passivity of the stack and of the film share over angle and thickness.
  for(double t:{0.1,1.0,2.0,3.0,6.0})for(double th=0;th<90;th+=2.5){Coeff f=film_stack(std::cos(th*PI/180),Z(2.2447212969283283,0.06982232081911283),t,Z(4.308219763636363,20.63042690909091),10.5);
    if(f.Ts< -1e-14||f.Tp< -1e-14||1-std::norm(f.rs)-f.Ts< -1e-12||1-std::norm(f.rp)-f.Tp< -1e-12)throw std::runtime_error("film stack not passive");}
  // (c) a lossless film (k=0) on Ti: all floor loss must be in the substrate.
  {Coeff f=film_stack(std::cos(40*PI/180),Z(1.5,0),2.0,Z(4.308219763636363,20.63042690909091),10.5);if(std::abs(1-std::norm(f.rs)-f.Ts)>1e-12||std::abs(1-std::norm(f.rp)-f.Tp)>1e-12)throw std::runtime_error("lossless film absorbs");}
  // (d) flat M3 face at normal incidence: tallies equal the analytic stack values.
  c.surface="M3";c.n=4.308219763636363;c.k=20.63042690909091;c.film_n=2.2447212969283283;c.film_k=0.06982232081911283;c.film_thickness=2;c.wavelength=10.5;
  {auto out3=trace(active,{0,0,0},d,initial,c);Coeff f=film_stack(1,Z(c.film_n,c.film_k),2,Z(c.n,c.k),10.5);
    for(int j=0;j<2;++j){double A=1-(j?std::norm(f.rp):std::norm(f.rs)),T=j?f.Tp:f.Ts;
      if(std::abs(out3.energy[j][ACTIVE]-A)>1e-13||std::abs(out3.energy[j][TIFLOOR]-T)>1e-13||std::abs(out3.energy[j][FILM]-(A-T))>1e-13)throw std::runtime_error("flat M3 tallies failed");}}
  // (e) M3 energy closure and sub-tally consistency in an A3 cavity at oblique incidence.
  {Geometry a3("A3",1,77,16.4,69.5455842631299,80);V dd=unit({-.35,-.2,-.91});auto init=initial_fields(std::atan2(.2,.35),dd);
    for(int j=0;j<300;++j){auto r=trace(a3,a3.sample(rng),dd,init,c);for(const auto& en:r.energy){
      if(std::abs(en[TOTAL]+en[ESCAPE]+en[RESIDUAL]-1)>1e-11)throw std::runtime_error("M3 energy closure failed");
      if(std::abs(en[FILM]+en[TIFLOOR]-en[ACTIVE])>1e-13||en[LAND]>en[WALL]+1e-15||en[FILM]< -1e-13||en[TIFLOOR]< -1e-13)throw std::runtime_error("M3 sub-tally inconsistency");}}}
  test_new_geometries();
  std::cout<<"{\"self_test\":\"PASS\",\"checks\":[\"normal_Fresnel\",\"normal_global_Jones_sign\",\"active_geometry_areas\",\"M0_M1_M2_energy_closure_900_rays\",\"grey_absorptance\",\"M3_zero_thickness_equals_bare_Fresnel_8_angles_3_films\",\"M3_stack_passivity_5_thicknesses_36_angles\",\"M3_lossless_film_no_film_loss\",\"M3_flat_normal_tallies\",\"M3_A3_energy_closure_and_subtallies_300_rays\",\"five_family_M0_MG_closure_and_direct_normal_catch\",\"mirrored_offset_and_uniform_scale_ray_pairs\",\"infinite_groove_y_transport\"],\"max_split_error\":"<<max_split_error<<"}\n";
}
} // namespace
int main(int argc,char**argv){try{
  Config c=args(argc,argv);std::cout<<std::setprecision(17);
  if(c.self_test){self_test();return 0;}if(c.fresnel_test){fresnel_output(c);return 0;}if(c.jones_probe){jones_output(c);return 0;}if(c.film_probe){film_output(c);return 0;}
  auto began=std::chrono::steady_clock::now();Geometry g(c.structure,c.segments,c.height,c.bottom_width,c.upper_width,c.pitch,c.offset_x);RNG rng{c.seed};double th=c.theta*PI/180,ph=c.phi*PI/180;
  V d={-std::sin(th)*std::cos(ph),-std::sin(th)*std::sin(ph),-std::cos(th)};auto initial=initial_fields(ph,d);
  std::array<Stats,3> stats;uint64_t bounces=0,wraps=0,cutoffs=0,caps=0,lost=0;int max_seen=0;double max_closure_error=0;
  for(uint64_t i=0;i<c.rays;++i){Result r=trace(g,g.sample(rng),d,initial,c);Energy avg{};
    for(int j=0;j<2;++j){stats[j].add(r.energy[j]);for(int k=0;k<METRICS;++k)avg[k]+=.5*r.energy[j][k];max_closure_error=std::max(max_closure_error,std::abs(r.energy[j][TOTAL]+r.energy[j][ESCAPE]+r.energy[j][RESIDUAL]-1));}
    stats[2].add(avg);bounces+=r.bounces;wraps+=r.wraps;max_seen=std::max(max_seen,r.bounces);cutoffs+=r.reason==1;caps+=r.reason==2;lost+=r.reason==3;
  }
  double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count();
  std::cout<<"{\"schema\":1,\"model\":\"OPAQUE_POLARIZED_GEOMETRIC_RAYS_NOT_FULL_WAVE\",\"surface\":\""<<c.surface<<"\",\"structure\":\""<<c.structure<<"\",\"theta_deg\":"<<c.theta<<",\"phi_deg\":"<<c.phi<<",\"wavelength_um\":"<<c.wavelength<<",\"n\":"<<c.n<<",\"k\":"<<c.k<<",\"grey_alpha\":"<<c.grey_alpha<<",\"cutoff\":"<<c.cutoff<<",\"rays\":"<<c.rays<<",\"seed\":"<<c.seed<<",\"curve_segments\":"<<c.segments<<",\"max_bounces\":"<<c.max_bounces<<",\"channels\":{\"s\":";
  print_stats(stats[0]);std::cout<<",\"p\":";print_stats(stats[1]);std::cout<<",\"unpolarized\":";print_stats(stats[2]);
  double active_area=g.active_area();
  std::cout<<"},\"ray_origin_policy\":\"exact_surface_barycentric_scale_roundoff_self_contact_v2\""
           <<",\"max_ray_closure_error\":"<<max_closure_error
           <<",\"cutoff_terminated_rays\":"<<cutoffs<<",\"cap_terminated_rays\":"<<caps<<",\"geometry_unresolved_rays\":"<<lost
           <<",\"mean_bounces\":"<<double(bounces)/c.rays<<",\"max_observed_bounces\":"<<max_seen
           <<",\"periodic_wraps\":"<<wraps<<",\"wall_seconds\":"<<seconds
           <<",\"unit_cell_area_um2\":"<<g.cell_area<<",\"active_area_um2\":"<<active_area
           <<",\"opening_area_um2\":"<<g.opening_area<<",\"triangles\":"<<g.mesh.tris.size()
           <<",\"height_um\":"<<c.height<<",\"bottom_width_um\":"<<c.bottom_width
           <<",\"upper_width_um\":"<<c.upper_width<<",\"pitch_um\":"<<c.pitch
           <<",\"offset_x_um\":"<<g.offset_x<<",\"width_definition\":\""<<g.width_definition()<<'"';
  if(c.surface=="M0"||c.surface=="MG")
    std::cout<<",\"floor_alpha\":"<<c.floor_alpha<<",\"wall_alpha\":"<<c.wall_alpha;
  else
    std::cout<<",\"floor_alpha\":null,\"wall_alpha\":null";
  std::cout<<",\"geometry_meta\":{\"profile\":\""<<g.profile()
           <<"\",\"opening_sides\":"<<g.opening_sides
           <<",\"footprint_type\":\""<<g.footprint_type()
           <<"\",\"floor_area_fraction\":"<<active_area/g.cell_area
           <<",\"opening_area_fraction\":"<<g.opening_area/g.cell_area
           <<",\"y_invariant\":"<<(g.trench?"true":"false")
           <<",\"mesh_validation\":\""<<g.mesh_validation()<<"\"}"
           <<",\"mesh_validation\":\""<<g.mesh_validation()<<'"'
           <<",\"curve_family\":\""<<(g.hex&&g.name!="A3"?"original_normalized_section_remapped_between_endpoints_not_all_ideal_CPC":g.profile())<<'"'
           <<",\"lattice\":\""<<(g.hex||g.round?"triangular_parameterized_primitive":(g.square||g.trench?"square_parameterized_primitive":"square_100um"))<<'"'
           <<",\"geometry_units\":\"um\",\"polarization\":\"paired_s_p_complex_global_Jones\""
           <<",\"grey_residual_phase\":\"PEC_like_rs_negative_rp_positive\""
           <<",\"uncertainty\":\"MC_SE_for_spatial_sampling_only_residual_is_separate_bound\""
           <<",\"film_n\":"<<c.film_n<<",\"film_k\":"<<c.film_k<<",\"film_thickness_um\":"<<c.film_thickness
           <<",\"max_layered_split_error\":"<<max_split_error
           <<",\"limitations\":\""
           <<(c.surface=="M3"?"M3_planar_coherent_single_layer_floor_on_semi_infinite_Ti_local_cosine_no_lateral_film_edges_no_buried_wall_no_diffraction_or_cavity_interference"
              :(c.surface=="MG"||c.surface=="M0"?"prescribed_opaque_specular_surface_no_spectrum_temperature_scattering_or_diffraction"
                :"opaque_local_bulk_Fresnel_no_diffraction_interference_roughness_film_or_temperature_model"))
           <<"\"}\n";
  return 0;
}catch(const std::exception&e){std::cerr<<"ray_material error: "<<e.what()<<'\n';return 2;}}
