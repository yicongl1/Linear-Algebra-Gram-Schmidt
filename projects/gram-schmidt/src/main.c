//  main.c
//  matrix
//  Created by 李毅聪 on 2022/2/11.

#include <stdio.h>
#include <math.h>
void norme1(double a[],int l)//该子程序用于将列向量转化为单位向量
{
    double s=0;
    for(int i=1;i<=l;i++)s+=a[i]*a[i];
    s=sqrt(s);
    for(int i=1;i<=l;i++)a[i]/=s;
}

double produitscalaire(double a[],double b[],int l)//该子程序用于计算两向量内积
{
    double s=0;
    for(int i=1;i<=l;i++)s+=a[i]*b[i];
    return s;
}

int main()
{
    int x,y;
    double a[100][100]={},b[100][100]={};//两个二维数组，前一个用于读入，后一个用于输出
    scanf("%d%d",&x,&y);
    for(int i=1;i<=y;i++)
        for(int j=1;j<=x;j++)
            scanf("%lf",a[j]+i);//转置输入，便于后面数组计算
    for(int i=1;i<=y;i++)
    {
        for(int j=1;j<=x;j++)b[i][j]=a[i][j];//完成uk=vk的第一步计算
        for(int j=1;j<i;j++)
        {
            double ps=produitscalaire(a[i],b[j],x);
            for(int k=1;k<=x;k++)
                b[i][k]-=ps*b[j][k];//每一项减去vk在uj上的投影
        }
        norme1(b[i],x);//得到同方向的单位向量
    }
    for(int i=1;i<=y;i++)
    {
        for(int j=1;j<=x;j++)
            printf("%10lf",b[j][i]);//转置输出，控制格式
        printf("\n");
    }
    return 0;
}
