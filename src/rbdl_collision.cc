// 本程序实现基于RBDL的刚体碰撞检测机制
#include "rbdl/rbdl_collision.h"

using namespace GeometryCollision;



double Clamp(double t, double A, double B)
{
    if (t < A)
    {
        return A;
    }
    else if (t > B)
    {
        return B;
    }
    else
    {
        return t;
    }
}

Vector3d ClosestPointOnLineSegment(Vector3d A, Vector3d B, Vector3d Point)
{
    Vector3d AB = B - A;
    //   float t = dot(Point – A, AB) / dot(AB, AB);
    Vector3d AP = Point - A;
    float t = AP.dot(AB);
    t = t / AB.norm();
    if (t < 0)
    {
        t = 0;
    }
    if (t > 1)
    {
        t = 1;
    }
    return A + t * AB; // saturate(t) can be written as: min((max(t, 0), 1)
}

// Computes closest points C1 and C2 of S1(s)=P1+s*(Q1-P1) and
// S2(t)=P2+t*(Q2-P2), returning s and t. Function result is squared
// distance between between S1(s) and S2(t)
double ClosestPtSegmentSegment(Vector3d p1, Vector3d q1, Vector3d p2, Vector3d q2,
                               double &s, double &t, Vector3d &c1, Vector3d &c2)
{
    Vector3d d1 = q1 - p1; // Direction vector of segment S1
    Vector3d d2 = q2 - p2; // Direction vector of segment S2
    Vector3d r = p1 - p2;
    double a = d1.dot(d1); // Squared length of segment S1, always nonnegative
    double e = d2.dot(d2); // Squared length of segment S2, always nonnegative
    double f = d2.dot(r);  //(d2, r);
    // Check if either or both segments degenerate into points
    if (a <= EPSILON && e <= EPSILON)
    {
        // Both segments degenerate into points
        s = t = 0.0f;
        c1 = p1;
        c2 = p2;
        return (c1 - c2).dot(c1 - c2); // Dot(c1 - c2, c1 - c2);
    }
    if (a <= EPSILON)
    {
        // First segment degenerates into a point
        s = 0.0f;
        t = f / e; // s = 0 => t = (b*s + f) / e = f / e
        t = Clamp(t, 0.0f, 1.0f);
    }
    else
    {
        double c = d1.dot(r); //(d1, r);
        if (e <= EPSILON)
        {
            // Second segment degenerates into a point
            t = 0.0;
            s = Clamp(-c / a, 0.0f, 1.0f); // t = 0 => s = (b*t - c) / a = -c / a
        }
        else
        {
            // The general nondegenerate case starts here
            double b = d1.dot(d2);        // Dot(d1, d2);
            double denom = a * e - b * b; // Always nonnegative
            // If segments not parallel, compute closest point on L1 to L2 and
            // clamp to segment S1. Else pick arbitrary s (here 0)
            if (denom != 0.0f)
            {
                s = Clamp((b * f - c * e) / denom, 0.0f, 1.0f);
            }
            else
                s = 0.0f;
            // Compute point on L2 closest to S1(s) using
            // t = Dot((P1 + D1*s) - P2,D2) / Dot(D2,D2) = (b*s + f) / e
            t = (b * s + f) / e;
            // If t in [0,1] done. Else clamp t, recompute s for the new value
            // of t using s = Dot((P2 + D2*t) - P1,D1) / Dot(D1,D1)= (t*b - c) / a
            // and clamp s to [0, 1]
            if (t < 0.0f)
            {
                t = 0.0f;
                s = Clamp(-c / a, 0.0f, 1.0f);
            }
            else if (t > 1.0f)
            {
                t = 1.0f;
                s = Clamp((b - c) / a, 0.0f, 1.0f);
            }
        }
    }
    c1 = p1 + d1 * s;
    c2 = p2 + d2 * t;
    return (c1 - c2).dot(c1 - c2); // Dot(c1 - c2, c1 - c2);
}

// Returns the squared distance between point c and segment ab
double SqDistPointSegment(Vector3d a, Vector3d b, Vector3d c)
{
    Vector3d ab = b - a;
    Vector3d ac = c - a;
    Vector3d bc = c - b;
    double e = ac.dot(ab); // Dot(ac, ab);
    // Handle cases where c projects outside ab
    if (e <= 0.0f)
        return ac.dot(ac); //;
    double f = ab.dot(ab); // Dot(ab, ab);
    if (e >= f)
        return bc.dot(bc); // Dot(bc, bc);
    // Handle cases where c projects onto ab
    return ac.dot(ac) - e * e / f; // Dot(ac, ac) – e * e / f;
}
double DistPointBox(Vector3d v, Vector3d sv)
// v和sv需要在第一像限,v表示box的尖，sv表示点
{
    double ret = 0;
    double temp1, temp2;
    if (sv[0] <= v[0] && sv[1] <= v[1] && sv[2] <= v[2]) // 圆心在box内部
    {
        ret = -(v - sv).norm();
    }
    else if (sv[0] > v[0] && sv[1] <= v[1] && sv[2] <= v[2]) // 圆心yz在内部，x在外部
    {
        ret = sv[0] - v[0];
    }
    else if (sv[0] <= v[0] && sv[1] > v[1] && sv[2] <= v[2])
    {
        ret = sv[1] - v[1];
    }
    else if (sv[0] <= v[0] && sv[1] <= v[1] && sv[2] > v[2])
    {
        ret = sv[2] - v[2];
    }
    else if (sv[0] > v[0] && sv[1] > v[1] && sv[2] <= v[2]) //
    {
        temp1 = sv[0] - v[0];
        temp1 = temp1 * temp1;
        temp2 = sv[1] - v[1];
        temp2 = temp2 * temp2;
        ret = sqrt(temp1 + temp2);
    }
    else if (sv[0] > v[0] && sv[1] <= v[1] && sv[2] > v[2]) //
    {
        temp1 = sv[0] - v[0];
        temp1 = temp1 * temp1;
        temp2 = sv[2] - v[2];
        temp2 = temp2 * temp2;
        ret = sqrt(temp1 + temp2);
    }
    else if (sv[0] <= v[0] && sv[1] > v[1] && sv[2] > v[2]) //
    {
        temp1 = sv[1] - v[1];
        temp1 = temp1 * temp1;
        temp2 = sv[2] - v[2];
        temp2 = temp2 * temp2;
        ret = sqrt(temp1 + temp2);
    }
    else
    {
        ret = (sv - v).norm();
    }
    return ret;
}
double DistPointCylinder(double cy_r,double cy_h, Vector3d sv)
{
    double ret = UNSUPPORT_RETURN;
    double p_rr=sv[0]*sv[0]+sv[1]*sv[1];//球体距离Z轴线距离平方
    double rr=cy_r*cy_r;//圆柱半径平方
    Vector3d edge_point(0,0,0);//圆柱上最近点
    if(sv[2]>=cy_h)//点在圆柱上方
    {
        if(p_rr<rr)
        {
            ret=sv[2]-cy_h;
        }
        else
        {
            edge_point[0]=sv[0]/sqrt(p_rr)*cy_r;
            edge_point[1]=sv[1]/sqrt(p_rr)*cy_r;
            edge_point[2]=cy_h;
            ret=(sv-edge_point).norm();
        }
    }
    else if(sv[2]<=0)//点在圆柱下方
    {
        if(p_rr<rr)
        {
            ret=-sv[2];
        }
        else
        {
            edge_point[0]=sv[0]/sqrt(p_rr)*cy_r;
            edge_point[1]=sv[1]/sqrt(p_rr)*cy_r;
            edge_point[2]=0;
            ret=(sv-edge_point).norm();
        }
    }
    else//点在圆柱腰部
    {
        if(p_rr<rr)//在圆柱内部
        {
            ret=min(min(sv[2]-cy_h,sqrt(p_rr)-cy_r),-sv[2]);
        }
        else//在圆柱外部
        {
            ret=sqrt(p_rr)-cy_r;
        }
    }
    return ret;
}


/******************************************************************************/

double Sphere_Sphere(Sphere *A, Sphere *B)
{
    // A.center_pos_base
    double ret = (A->pos_Base - B->pos_Base).norm() - A->radius - B->radius;
    return ret;
}
double Sphere_Box(Sphere *A, Box *B)
{
    // A.center_pos_base
    // A->T_Base_Geo.r-B->T_Base_Geo.r

    // A->T_Base_Geo.r
    // B->T_Base_Geo.E将Base的向量转到
    double ret = 0;
    Vector3d sphere_loc_in_B = B->T_Base_Geo.E * A->pos_Base - B->T_Base_Geo.E * (B->T_Base_Geo.r); // 计算球圆心在Box坐标系中的表示
    Vector3d sv;
    sv[0] = abs(sphere_loc_in_B[0]);
    sv[1] = abs(sphere_loc_in_B[1]);
    sv[2] = abs(sphere_loc_in_B[2]);
    ret = DistPointBox(B->v, sv) - A->radius;
}
double Sphere_Capsule(Sphere *A, Capsule *B)
{
    double ret = 0;
    ret = SqDistPointSegment(B->pA_Base, B->pB_Base, A->pos_Base);
    return sqrt(ret) - A->radius - B->radius;
}
double Sphere_Plane(Sphere *A, Plane *B)
{
    printf("Sphere_Plane unsupport!");return UNSUPPORT_RETURN;
}
double Box_Box(Box *A, Box *B)
{
    printf("Box_Box unsupport!");return UNSUPPORT_RETURN;
}
double Box_Plane(Box *A, Plane *B)
{
    printf("Box_Plane unsupport!");return UNSUPPORT_RETURN;
}
double Box_Capsule(Box *A, Capsule *B)
{
    double ret = 1e20, ret_temp;
    Vector3d box_r = A->T_Base_Geo.E * (A->T_Base_Geo.r);
    Vector3d p1 = A->T_Base_Geo.E * B->pA_Base - box_r; // 计算capsule球圆心在Box坐标系中的表示
    Vector3d p2 = A->T_Base_Geo.E * B->pB_Base - box_r; // 计算capsule球圆心在Box坐标系中的表示
    Vector3d sv;

    double len = (p1 - p2).norm();
    Vector3d p12 = (p2 - p1) / len;

    // double t=0;
    // double delt_t=floor(len/(B->radius*B->capsule_radius_div));
    int i = 0;
    double delt_l = B->radius * B->capsule_radius_div;
    int imax = floor(len / delt_l);
    do
    {
        sv = p1 + p12 * delt_l * i;
        sv[0] = abs(sv[0]);
        sv[1] = abs(sv[1]);
        sv[2] = abs(sv[2]);
        ret_temp = DistPointBox(A->v, sv) - B->radius;
        if (ret_temp < ret)
        {
            ret = ret_temp;
        }
        i++;
    }while (i < imax);

    p2[0] = abs(p2[0]);
    p2[1] = abs(p2[1]);
    p2[2] = abs(p2[2]);
    ret_temp = DistPointBox(A->v, p2) - B->radius;
    if (ret_temp < ret)
    {
        ret = ret_temp;
    }
    return ret;

    
}
double Capsule_Capsule(Capsule *A, Capsule *B)
{
    double ret = 0;
    double s, t;
    Vector3d A3, B3;
    ret = ClosestPtSegmentSegment(A->pA_Base, A->pB_Base, B->pA_Base, B->pB_Base, s, t, A3, B3);
    ret = sqrt(ret) - A->radius - B->radius;
    return ret;
}
double Capsule_Plane(Capsule *A, Plane *B)
{
    printf("Capsule_Plane unsupport!");return UNSUPPORT_RETURN;
}
double Plane_Plane(Plane *A, Plane *B)
{
    printf("Plane_Plane unsupport!");return UNSUPPORT_RETURN;
}
double Cylinder_Sphere(Sphere *A,Cylinder *B)
{
    double ret = UNSUPPORT_RETURN;
    Vector3d sphere_loc_in_B = B->T_Base_Geo.E * A->pos_Base - B->T_Base_Geo.E * (B->T_Base_Geo.r); // 计算球圆心在Cylinder坐标系中的表示
    ret=DistPointCylinder(B->radius,B->height, sphere_loc_in_B)-A->radius;
    return ret;
}
double Cylinder_Box(Box *A, Cylinder *B)
{
    printf("Cylinder_Box unsupport!");return UNSUPPORT_RETURN;
}
double Cylinder_Capsule(Capsule *A,Cylinder *B)
{
    printf("Cylinder_Capsule unsupport!");return UNSUPPORT_RETURN;
}
double Cylinder_Plane(Plane *A,Cylinder *B)
{
    printf("Cylinder_Plane unsupport!");return UNSUPPORT_RETURN;
}
double Cylinder_Cylinder(Cylinder *A, Cylinder *B)
{
    printf("Cylinder_Cylinder unsupport!");return UNSUPPORT_RETURN;
}

/******************************************************************************/

Sphere::Sphere(Vector3d pos_, double radius_)
{
    radius = radius_;
    pos = pos_;
    geometry_type = GEOMETRY_SPHERE;
};
double Sphere::collisionDistance(GeometryBase *geo)
{
    double min_dis = 0;
    switch (geo->geometry_type)
    {
    case GEOMETRY_SPHERE: // 圆与圆
        min_dis = Sphere_Sphere(this, (Sphere *)(geo));
        break;
    case GEOMETRY_BOX:
        min_dis = Sphere_Box(this, (Box *)(geo));
        break;
    case GEOMETRY_CAPSULE:
        min_dis = Sphere_Capsule(this, (Capsule *)(geo));
        break;
    case GEOMETRY_PLANE:
        min_dis = Sphere_Plane(this, (Plane *)(geo));
        break;
    case GEOMETRY_CYLINDER:
        min_dis = Cylinder_Sphere(this, (Cylinder *)(geo));
        break;
    default:
        break;
    }
    return min_dis;
};
void Sphere::print_data()
{
    printf("Sphere:pos(%f,%f,%f),rad(%f)##", this->pos_Base[0], this->pos_Base[1], this->pos_Base[2], this->radius);
}
void Sphere::update_data(Math::SpatialTransform T_Base_body)
{
    pos_Base = (T_Base_body.E.transpose()) * pos + T_Base_body.r;
}

Box::Box(Math::SpatialTransform T_Body_Geo_, double x, double y, double z)
{
    v[0] = x / 2;
    v[1] = y / 2;
    v[2] = z / 2;
    T_Body_Geo = T_Body_Geo_;
    geometry_type = GEOMETRY_BOX;
};
double Box::collisionDistance(GeometryBase *geo)
{
    double min_dis = 0;
    switch (geo->geometry_type)
    {
    case GEOMETRY_SPHERE: // Box与圆
        min_dis = Sphere_Box((Sphere *)(geo), this);
        break;
    case GEOMETRY_BOX:
        min_dis = Box_Box((Box *)(geo), this);
        break;
    case GEOMETRY_CAPSULE:
        min_dis = Box_Capsule(this, (Capsule *)(geo));
        break;
    case GEOMETRY_PLANE:
        min_dis = Box_Plane(this, (Plane *)(geo));
        break;
    case GEOMETRY_CYLINDER:
        min_dis = Cylinder_Box(this, (Cylinder *)(geo));
        break;
    default:
        break;
    }
    return min_dis;
};
void Box::print_data()
{
    printf("Box:pos(%f,%f,%f),xyz(%f,%f,%f)##", this->T_Base_Geo.r[0], this->T_Base_Geo.r[1], this->T_Base_Geo.r[2], this->v[0] * 2, this->v[1] * 2, this->v[2] * 2);
}
void Box::update_data(Math::SpatialTransform T_Base_body)
{
    // geometry_list[i]->T_Base_Geo=geometry_list[i]->T_Body_Geo*T_Base_Body;
    T_Base_Geo = T_Body_Geo * T_Base_body;
}

Capsule::Capsule(Vector3d pA_, Vector3d pB_, double radius_)
{
    pA = pA_;
    pB = pB_;
    z_half = (pA_ - pB_).norm() / 2;
    radius = radius_;
    geometry_type = GEOMETRY_CAPSULE;
}
double Capsule::collisionDistance(GeometryBase *geo)
{
    double min_dis = 0;
    switch (geo->geometry_type)
    {
    case GEOMETRY_SPHERE:
        min_dis = Sphere_Capsule((Sphere *)geo, this);
        break;
    case GEOMETRY_BOX:
        min_dis = Box_Capsule((Box *)geo, this);
        break;
    case GEOMETRY_CAPSULE:
        min_dis = Capsule_Capsule((Capsule *)geo, this);
        break;
    case GEOMETRY_PLANE:
        min_dis = Capsule_Plane(this, (Plane *)(geo));
        break;
    case GEOMETRY_CYLINDER:
        min_dis = Cylinder_Capsule(this, (Cylinder *)(geo));
        break;
    default:
        break;
    }
    return min_dis;
};
void Capsule::print_data()
{
    printf("Capsule:Apos(%f,%f,%f),Bpos(%f,%f,%f),rad(%f)##",
           pA_Base[0], pA_Base[1], pA_Base[2],
           pB_Base[0], pB_Base[1], pB_Base[2],
           this->radius);
}
void Capsule::update_data(Math::SpatialTransform T_Base_body)
{
    Math::Matrix3d E_ = T_Base_body.E.transpose();
    pA_Base = E_ * pA + T_Base_body.r;
    pB_Base = E_ * pB + T_Base_body.r;
}

Plane::Plane(Vector3d p_, Vector3d n_)
{
    p = p_;
    n = n_;
    geometry_type = GEOMETRY_PLANE;
}
double Plane::collisionDistance(GeometryBase *geo)
{
    double min_dis = 0;
    switch (geo->geometry_type)
    {
    case GEOMETRY_SPHERE:
        min_dis = Sphere_Plane((Sphere *)geo, this);
        break;
    case GEOMETRY_BOX:
        min_dis = Box_Plane((Box *)geo, this);
        break;
    case GEOMETRY_CAPSULE:
        min_dis = Capsule_Plane((Capsule *)geo, this);
        break;
    case GEOMETRY_PLANE:
        min_dis = Plane_Plane((Plane *)geo, this);
        break;
    case GEOMETRY_CYLINDER:
        min_dis = Cylinder_Plane(this, (Cylinder *)(geo));
        break;
    default:
        break;
    }
    return min_dis;
};
void Plane::print_data()
{
    printf("Plane##");
}
void Plane::update_data(Math::SpatialTransform T_Base_body)
{
    return;
}

Cylinder::Cylinder(Math::SpatialTransform T_Body_Geo_,double r_,double h_)
{
    T_Body_Geo=T_Body_Geo_;
    radius=r_;
    height=h_;
}
void Cylinder::print_data()
{
    printf("Cylinder:Apos(%f,%f,%f),Bpos(%f,%f,%f),rad(%f)##",
        pA_Base[0], pA_Base[1], pA_Base[2],
        pB_Base[0], pB_Base[1], pB_Base[2],
        this->radius);
}
void Cylinder::update_data(Math::SpatialTransform T_Base_body)
{
    T_Base_Geo = T_Body_Geo * T_Base_body;
    Math::Matrix3d E_ = T_Base_Geo.E.transpose();
    pA_Base = T_Base_Geo.r;
    v_Base=E_*Math::Vector3d(0,0,1);
    pB_Base = T_Base_Geo.r+v_Base*height;

    
}
double Cylinder::collisionDistance(GeometryBase *geo)
{
    double min_dis = 0;
    switch (geo->geometry_type)
    {
    case GEOMETRY_SPHERE: // Box与圆
        min_dis = Cylinder_Sphere((Sphere *)(geo), this);
        break;
    case GEOMETRY_BOX:
        min_dis = Cylinder_Box((Box *)(geo), this);
        break;
    case GEOMETRY_CAPSULE:
        min_dis = Cylinder_Capsule((Capsule *)(geo),this);
        break;
    case GEOMETRY_PLANE:
        min_dis = Cylinder_Plane((Plane *)(geo),this);
        break;
    case GEOMETRY_CYLINDER:
        min_dis = Cylinder_Cylinder(this, (Cylinder *)(geo));
        break;
    default:
        break;
    }
    return min_dis;
};

/*****************************************/

GeometryBase::GeometryBase()
{
    geometry_type = GEOMETRY_BASE;
}
double GeometryBase::collisionDistance(GeometryBase *geo)
{
    printf("collision check unsupport!");return UNSUPPORT_RETURN;
}
void GeometryBase::print_data()
{
    printf("GeometryBase ");
}
void GeometryBase::update_data(Math::SpatialTransform T_Base_body)
{
    return;
}

/***********************************************************/

void GeometryGroup::updatGeometries()
{
    for (int i = 0; i < geometry_list.size(); i++)
    {
        geometry_list[i]->update_data(T_Base_Body);
    }
};

double GeometryGroup::collisionDistanceCalc(GeometryGroup &gp)
{
    double min_dis = 1e20, dis_temp;
    int i = 0, j = 0;
    int si, sj;
    si = geometry_list.size();
    sj = gp.geometry_list.size();
    for (i = 0; i < si; i++)
    {
        for (j = 0; j < sj; j++)
        {
            dis_temp = geometry_list[i]->collisionDistance(gp.geometry_list[j]);
            if (dis_temp < min_dis)
            {
                min_dis = dis_temp;
            }
        }
    }
    return min_dis;
}

double GeometryGroup::collisionDistanceCalc(GeometryBase &geo)
{
    double min_dis = 1e20, dis_temp;
    int i = 0, j = 0;
    int si, sj;
    si = geometry_list.size();

    for (i = 0; i < si; i++)
    {
        dis_temp = geometry_list[i]->collisionDistance(&geo);
        if (dis_temp < min_dis)
        {
            min_dis = dis_temp;
        }
    }
    return min_dis;
}