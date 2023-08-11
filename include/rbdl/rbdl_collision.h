//本程序实现基于RBDL的刚体碰撞检测机制
#include <iostream>
#include <limits>
#include <cstring>
#include <assert.h>

#include "rbdl/rbdl_mathutils.h"
#include "SpatialAlgebraOperators.h"
#include <iostream> 
#include <algorithm> 
#include <vector> 
using namespace std; 
//RBDL_DLLAPI void 
using namespace RigidBodyDynamics;
using namespace RigidBodyDynamics::Math;



namespace GeometryCollision
{

#define UNSUPPORT_RETURN 1e20 //如果不支持碰撞检测，则返回的默认值
#define EPSILON 1e-8
#define VERY_LARGE UNSUPPORT_RETURN

enum GeometryTypes
{
    GEOMETRY_BASE,
    GEOMETRY_SPHERE,
    GEOMETRY_BOX,
    GEOMETRY_CAPSULE,
    GEOMETRY_PLANE,
    GEOMETRY_CYLINDER
};

class GeometryBase//几何体的基类
{
    public:
    GeometryBase();
    //Math::SpatialTransform T_Body_Geo;//transform from body frame to geometry frame
    //Math::SpatialTransform T_Base_Geo;//transform from base frame to geometry frame
    int geometry_type;
    virtual void print_data();
    virtual double collisionDistance(GeometryBase *geo);
    virtual void update_data(Math::SpatialTransform T_Base_body);
};

class Sphere:public GeometryBase//任意位置球体
{
    public:
    Sphere(Vector3d pos_,double radius_);
    double collisionDistance(GeometryBase *geo);
    void print_data();
    double radius;//半径
    Vector3d pos;//在物体坐标的位置
    Vector3d pos_Base;//在基坐标的位置
    void update_data(Math::SpatialTransform T_Base_body);
    
};
class Box:public GeometryBase//任意姿态位置的长方体
{
    public:
    Box(Math::SpatialTransform T_Body_Geo_,double x,double y,double z);//位姿，长，宽，高
    double collisionDistance(GeometryBase *geo);
    void print_data();
    Math::SpatialTransform T_Body_Geo;
    Math::SpatialTransform T_Base_Geo;
    Vector3d v;
    void update_data(Math::SpatialTransform T_Base_body);
};
class Capsule:public GeometryBase//中心为AB的胶囊体
{
    public:
    Capsule(Vector3d pA_,Vector3d pB_,double radius_);
    double collisionDistance(GeometryBase *geo);
    void print_data();
    Vector3d pA,pA_Base;
    Vector3d pB,pB_Base;
    double z_half;//
    double radius;
    double capsule_radius_div=1; //计算capsule和box的碰撞时，用多个球代替cpsule,capsule_radius_div表示每隔capsule_radius_div*radius距离就建立个点圆
    void update_data(Math::SpatialTransform T_Base_body);
    
};
class Cylinder:public GeometryBase//底部与Geo坐标系原点重合，方向沿Z轴，半径r，高度h
{
    public:
    Cylinder(Math::SpatialTransform T_Body_Geo_,double r_,double h_);//方向沿Z轴，半径r，高度h
    double collisionDistance(GeometryBase *geo);
    void print_data();
    Math::SpatialTransform T_Body_Geo;
    Math::SpatialTransform T_Base_Geo;
    double radius,height;
    Vector3d pA_Base;//cylinder的点
    Vector3d pB_Base;
    Vector3d v_Base;//cylinder在base坐标系的方向向量
    void update_data(Math::SpatialTransform T_Base_body);
};

class Plane:public GeometryBase//XY平面
{
    public:
    Plane(Vector3d p_,Vector3d n_);
    double collisionDistance(GeometryBase *geo);
    void print_data();
    Vector3d p,p_Base;//平面上的任意点
    Vector3d n,n_Base;//平面的法向量
    void update_data(Math::SpatialTransform T_Base_body);
};

class GeometryGroup //多个几何体的集合组成body
{
    public:
    Math::SpatialTransform T_Base_Body;//transform frame from base frame to group(body)frame
    std::vector <GeometryBase *> geometry_list;
    void updatGeometries();
    double collisionDistanceCalc(GeometryGroup &gp);//计算和另外外的几何grop的最小距离
    double collisionDistanceCalc(GeometryBase &geo);////计算和另外外的几何体的最小距离
};

};
