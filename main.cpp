#include<iostream>
#include<vector>
#include<random>
#include<glm/glm.hpp>
#include"svpng.inc"
#include<omp.h>	//openmp多线程加速

using namespace glm;
//
const int SAMPLE = 4096;
// 每次采样的亮度
const double BRIGHTNESS = (2.0f * 3.1415926f) * (1.0f / double(SAMPLE));

//输出图像分辨率
const int WIDTH = 256;
const int HEIGHT = 256;
//相机参数
const double SCREEN_Z = 1.1;
const glm::vec3 EYE = glm::vec3(0, 0, 4.0);   // 相机位置

void imShow(double* SRC);

// 颜色
const glm::vec3 RED(1, 0.5, 0.5);
const glm::vec3 GREEN(0.5, 1, 0.5);
const glm::vec3 BLUE(0.5, 0.5, 1);
const glm::vec3 YELLOW(1.0, 1.0, 0.1);
const glm::vec3 CYAN(0.1, 1.0, 1.0);
const glm::vec3 MAGENTA(1.0, 0.1, 1.0);
const glm::vec3 GRAY(0.5, 0.5, 0.5);
const glm::vec3 WHITE(1, 1, 1);
//光线
struct Ray
{
	glm::vec3 startPoint = glm::vec3(0,0,0);
	glm::vec3 direction = glm::vec3(0,0,0);
};

//物体表面材质定义
struct Material
{
	bool isEmissive = false;					//是否发光
	glm::vec3 normal = glm::vec3(0, 0, 0);		//法线
	glm::vec3 color = glm::vec3(0, 0, 0);		//颜色
	double specularRate = 0.0f;					//反射光占比
	double roughness = 1.0f;					//粗糙程度

	double refractRate = 0.0f;      // 折射光占比
	double refractAngle = 1.0f;     // 折射率
	double refractRoughness = 0.0f; // 折射粗糙度
};

//求交结果
struct HitResult
{
	bool isHit = false;
	double distance = 0.0;
	glm::vec3 hitPoint = glm::vec3(0, 0, 0);
	Material material;
};

//0-1随机数生成
double randf()
{
	static std::random_device rd;
	static std::mt19937 generator(rd());
	static std::uniform_real_distribution<double> distribution(0.0, 1.0);
	return distribution(generator);
}

//单位球内的随机向量
glm::vec3 randomVec3()
{
	glm::vec3 d;
	do
	{
		d = 2.0f * glm::vec3(randf(), randf(), randf()) - glm::vec3(1.f, 1.f, 1.f);	//(0~1)  -> (-1~1)
	} while (glm::dot(d, d) > 1.0f);
	return glm::normalize(d);
}

//法向半球随机向量
glm::vec3 randomDirection(glm::vec3 n)
{
	// 法向半球
	//glm::vec3 d;
	//do
	//{
	//	d = randomVec3();
	//} while (glm::dot(d, n) < 0.0f);
	//return d;
	return glm::normalize(randomVec3() + n);
}

class Shape
{
public:
	Shape() {}
	virtual HitResult intersect(Ray ray) { return HitResult(); }
};

//Triangle
class Triangle:public Shape
{
public:

	Triangle(){}
	Triangle(glm::vec3 P1, glm::vec3 P2, glm::vec3 P3, glm::vec3 color)
	{
		p1 = P1, p2 = P2, p3 = P3;
		material.normal = glm::normalize(glm::cross(p2 - p1, p3 - p1));
		material.color = color;
	}
	//与光线求交
	 HitResult intersect(Ray ray)
	{
		HitResult res;

		vec3 S = ray.startPoint;        // 射线起点
		vec3 d = ray.direction;         // 射线方向
		vec3 N = material.normal;       // 法向量
		if (dot(N, d) > 0.0f) N = -N;   // 获取正确的法向量

		// 如果视线和三角形平行
		if (fabs(dot(N, d)) < 0.00001f) return res;

		// 距离
		float t = (dot(N, p1) - dot(S, N)) / dot(d, N);
		if (t < 0.0005f) return res;    // 如果三角形在相机背面

		// 交点计算
		vec3 P = S + d * t;

		// 判断交点是否在三角形中
		vec3 c1 = cross(p2 - p1, P - p1);
		vec3 c2 = cross(p3 - p2, P - p2);
		vec3 c3 = cross(p1 - p3, P - p3);
		vec3 n = material.normal;   // 需要 "原生法向量" 来判断
		if (dot(c1, n) < 0 || dot(c2, n) < 0 || dot(c3, n) < 0) return res;

		// 装填返回结果
		res.isHit = true;
		res.distance = t;
		res.hitPoint = P;
		res.material = material;
		res.material.normal = N;    // 要返回正确的法向
		return res;
	}
public:
	glm::vec3 p1, p2, p3;
	Material material;
};

class Sphere :public Shape
{
public:

	Sphere() {}
	Sphere(vec3 o, double r, vec3 c) { O = o; R = r; material.color = c; }
	vec3 O;             // 圆心
	double R;           // 半径
	Material material;  // 材质
	// 与光线求交
	 HitResult intersect(Ray ray)
	{
		HitResult res;

		vec3 S = ray.startPoint;        // 射线起点
		vec3 d = ray.direction;         // 射线方向

		float OS = length(O - S);
		float SH = dot(O - S, d);
		float OH = sqrt(pow(OS, 2) - pow(SH, 2));

		if (OH > R) return res; // OH大于半径则不相交

		float PH = sqrt(pow(R, 2) - pow(OH, 2));

		float t1 = length(SH) - PH;
		float t2 = length(SH) + PH;
		float t = (t1 < 0) ? (t2) : (t1);   // 最近距离
		vec3 P = S + t * d;     // 交点

		// 防止自己交自己
		if (fabs(t1) < 0.0005f || fabs(t2) < 0.0005f) return res;

		// 装填返回结果
		res.isHit = true;
		res.distance = t;
		res.hitPoint = P;
		res.material = material;
		res.material.normal = normalize(P - O); // 要返回正确的法向
		return res;
	}
};



//输出SRC数组中的数据到图像
void imShow(double* SRC)
{
	unsigned char* image = new unsigned char[WIDTH * HEIGHT * 3];//image buffer
	unsigned char* p = image;
	double* srcImg = SRC;

	FILE* fp;
	fopen_s(&fp, "prt.png", "wb");

	for (int i = 0; i < HEIGHT; i++)
	{
		for (int j = 0; j < WIDTH; j++)
		{
			*p++ = (unsigned char)glm::clamp(pow(*srcImg++,1.0f/2.2f) * 255, 0.0, 255.0);	//R 通道
			*p++ = (unsigned char)glm::clamp(pow(*srcImg++,1.0f/2.2f) * 255, 0.0, 255.0);	//G 通道
			*p++ = (unsigned char)glm::clamp(pow(*srcImg++,1.0f/2.2f) * 255, 0.0, 255.0);	//B 通道
		}
	}
	//void svpng(FILE * fp, unsigned w, unsigned h, const unsigned char* img, int alpha)
	svpng(fp, WIDTH, HEIGHT, image, 0);
}

HitResult shoot(std::vector<Shape*> shapes, Ray ray)
{
	HitResult res, r;
	res.distance = 1145141919.810f;

	for (auto& shape : shapes)
	{
		r = shape->intersect(ray);
		if (r.isHit && r.distance < res.distance)	//记录距离最近的求交结果
			res = r;
	}
	return res;
}

//路径追踪
glm::vec3 pathTracing(std::vector<Shape*> shapes, Ray ray,int depth)
{
	if (depth > 8) return vec3(0);
	HitResult res = shoot(shapes, ray);

	if (!res.isHit) return vec3(0); // 未命中

	// 如果发光则返回颜色
	if (res.material.isEmissive) return res.material.color;

	// 有 P 的概率终止
	double r = randf();
	float P = 0.8;
	if (r > P) return vec3(0);

	// 否则继续
	Ray randomRay;
	randomRay.startPoint = res.hitPoint;
	randomRay.direction = randomDirection(res.material.normal);

	vec3 color = vec3(0);
	float cosine = fabs(dot(-ray.direction, res.material.normal));

	// 根据反射率决定光线最终的方向
	r = randf();
	if (r < res.material.specularRate)  // 镜面反射
	{
		vec3 ref = normalize(reflect(ray.direction, res.material.normal));
		randomRay.direction = mix(ref, randomRay.direction, res.material.roughness);
		color = pathTracing(shapes, randomRay, depth + 1) * cosine;
	}
	else if (res.material.specularRate <= r && r <= res.material.refractRate)    // 折射
	{
		vec3 ref = normalize(refract(ray.direction, res.material.normal, float(res.material.refractAngle)));
		randomRay.direction = mix(ref, -randomRay.direction, res.material.refractRoughness);
		color = pathTracing(shapes, randomRay, depth + 1) * cosine;
	}
	else    // 漫反射
	{
		vec3 srcColor = res.material.color;
		vec3 ptColor = pathTracing(shapes, randomRay, depth + 1) * cosine;
		color = ptColor * srcColor;    // 和原颜色混合
	}

	return color / P;
}


int main()
{
	std::vector<Shape*> shapes;
	Sphere s1 = Sphere(vec3(-0.65, -0.7, 0.0), 0.3, GREEN);
	Sphere s2 = Sphere(vec3(0.0, -0.3, 0.0), 0.4, WHITE);
	Sphere s3 = Sphere(vec3(0.65, 0.1, 0.0), 0.3, BLUE);
	s1.material.specularRate = 0.3;
	s1.material.roughness = 0.1;

	s2.material.specularRate = 0.3;
	s2.material.refractRate = 0.95;
	s2.material.refractAngle = 0.1;
	//s2.material.refractRoughness = 0.05;

	s3.material.specularRate = 0.3;

	shapes.push_back(&s1);
	shapes.push_back(&s2);
	shapes.push_back(&s3);

	shapes.push_back(new Triangle(vec3(-0.15, 0.4, -0.6), vec3(-0.15, -0.95, -0.6), vec3(0.15, 0.4, -0.6), YELLOW));
	shapes.push_back(new Triangle(vec3(0.15, 0.4, -0.6), vec3(-0.15, -0.95, -0.6), vec3(0.15, -0.95, -0.6), YELLOW));

	Triangle tt = Triangle(vec3(-0.2, -0.2, -0.95), vec3(0.2, -0.2, -0.95), vec3(-0.0, -0.9, 0.4), YELLOW);
	//tt.material.specularRate = 0.1;
	//tt.material.refractRate = 0.85;
	//tt.material.refractRoughness = 0.3;
	//shapes.push_back(&tt);

	// 发光物
	Triangle l1 = Triangle(vec3(0.4, 0.99, 0.4), vec3(-0.4, 0.99, -0.4), vec3(-0.4, 0.99, 0.4), WHITE);
	Triangle l2 = Triangle(vec3(0.4, 0.99, 0.4), vec3(0.4, 0.99, -0.4), vec3(-0.4, 0.99, -0.4), WHITE);
	l1.material.isEmissive = true;
	l2.material.isEmissive = true;
	shapes.push_back(&l1);
	shapes.push_back(&l2);

	// 背景盒子
	// bottom
	shapes.push_back(new Triangle(vec3(1, -1, 1), vec3(-1, -1, -1), vec3(-1, -1, 1), WHITE));
	shapes.push_back(new Triangle(vec3(1, -1, 1), vec3(1, -1, -1), vec3(-1, -1, -1), WHITE));
	// top
	shapes.push_back(new Triangle(vec3(1, 1, 1), vec3(-1, 1, 1), vec3(-1, 1, -1), WHITE));
	shapes.push_back(new Triangle(vec3(1, 1, 1), vec3(-1, 1, -1), vec3(1, 1, -1), WHITE));
	// back
	shapes.push_back(new Triangle(vec3(1, -1, -1), vec3(-1, 1, -1), vec3(-1, -1, -1), CYAN));
	shapes.push_back(new Triangle(vec3(1, -1, -1), vec3(1, 1, -1), vec3(-1, 1, -1), CYAN));
	// left
	shapes.push_back(new Triangle(vec3(-1, -1, -1), vec3(-1, 1, 1), vec3(-1, -1, 1), BLUE));
	shapes.push_back(new Triangle(vec3(-1, -1, -1), vec3(-1, 1, -1), vec3(-1, 1, 1), BLUE));
	// right
	shapes.push_back(new Triangle(vec3(1, 1, 1), vec3(1, -1, -1), vec3(1, -1, 1), RED));
	shapes.push_back(new Triangle(vec3(1, -1, -1), vec3(1, 1, 1), vec3(1, 1, -1), RED));


	double* image = new double[WIDTH * HEIGHT * 3];
	memset(image, 0.0, sizeof(double) * WIDTH * HEIGHT * 3);

	omp_set_num_threads(50); // 线程个数
	#pragma omp parallel for
	for (int k = 0; k < SAMPLE; k++)
	{
		double* p = image;

		for (int i = 0; i < HEIGHT; i++)
		{
			for (int j = 0; j < WIDTH; j++)
			{
				//像素坐标转投影平面坐标
				double x = 2.0 * double(j) / double(WIDTH) - 1.0;
				double y = 2.0 * double(HEIGHT - i) / double(HEIGHT) - 1.0;

				// MSAA
				x += (randf() - 0.5f) / double(WIDTH);
				y += (randf() - 0.5f) / double(HEIGHT);

				glm::vec3 coord = glm::vec3(x, y, SCREEN_Z);	//计算投影平面坐标
				glm::vec3 direction = glm::normalize(coord - EYE);	//计算光线投射方向

				// 生成光线
				Ray ray;
				ray.startPoint = coord;
				ray.direction = direction;

				//与场景的交点
				HitResult res = shoot(shapes, ray);
				glm::vec3 color = glm::vec3(0.0,0.0,0.0);

				if (res.isHit)
				{
					//命中光源直接返回颜色
					if (res.material.isEmissive)
						color = res.material.color;
					//命中实体则选择一个随机方向重新发射光线并进行路线追踪
					else
					{
						//根据交点处的法向量生成交点处反射的随机半球向量
						Ray randomRay; 
						randomRay.startPoint = res.hitPoint;
						randomRay.direction = randomDirection(res.material.normal);

						//根据反射率决定光线最终的方向

						double r = randf();
						if (r < res.material.specularRate)  // 镜面反射
						{
							vec3 ref = normalize(reflect(ray.direction, res.material.normal));
							randomRay.direction = mix(ref, randomRay.direction, res.material.roughness);
							color = pathTracing(shapes, randomRay, 0);
						}
						else if (res.material.specularRate <= r && r <= res.material.refractRate)    // 折射
						{
							vec3 ref = normalize(refract(ray.direction, res.material.normal, float(res.material.refractAngle)));
							randomRay.direction = mix(ref, -randomRay.direction, res.material.refractRoughness);
							color = pathTracing(shapes, randomRay, 0);
						}
						else    // 漫反射
						{
							vec3 srcColor = res.material.color;
							vec3 ptColor = pathTracing(shapes, randomRay, 0);
							color = ptColor * srcColor;    // 和原颜色混合
						}
						color *= BRIGHTNESS;
					}
				}
				*p += color.x; p++;
				*p += color.y; p++;
				*p += color.z; p++;
			}
		}
	}
	imShow(image);
	std::cout << "image output done.\n";
}