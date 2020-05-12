clear all;
clc;
%load data
Node=load('node_data.txt');
Line=load('impedance_map.txt');
n=size(Node,1);
%build Y
Y=zeros(n);
for ii=1:size(Line,1)
    Y(Line(ii,1),Line(ii,2))=-1/(Line(ii,3)+1i*Line(ii,4));
    Y(Line(ii,2),Line(ii,1))=-1/(Line(ii,3)+1i*Line(ii,4));
    Y(Line(ii,1),Line(ii,1))=Y(Line(ii,1),Line(ii,1))+1/(Line(ii,3)+1i*Line(ii,4))+1i*Line(ii,5);
    Y(Line(ii,2),Line(ii,2))=Y(Line(ii,2),Line(ii,2))+1/(Line(ii,3)+1i*Line(ii,4))+1i*Line(ii,5);
end

%initialize e f x
G=real(Y);
B=imag(Y);
e=zeros(n,1);
f=zeros(n,1);
x=zeros(2*n,1);
for a=1:n
    if Node(a,2)==1 %slack
        e(a)=1.1;
        f(a)=0;
        x(2*a-1)=e(a);
        x(2*a)=f(a);
    else %PVorPQ
        e(a)=1;
        f(a)=0;
        x(2*a-1)=e(a);
        x(2*a)=f(a);
    end
end

%initialize ys
ys=zeros(2*n,1);
for a=1:n
    if Node(a,2)==1 %slack
        ys(2*a-1)=0;
        ys(2*a)=0;
    end
    if Node(a,2)==2 %PV
        ys(2*a-1)=Node(a,3);
        ys(2*a)=Node(a,5)^2;
    end
    if Node(a,2)==3 %PQ
        ys(2*a-1)=-1*Node(a,3);
        ys(2*a)=-1*Node(a,4);
    end
end
  
maxtime=20;
for time=1:maxtime
    times=['times=' num2str(time)];
    disp(times);
    %y1 y2 y
    y1=zeros(n,1);
    y2=zeros(n,1);
    y=zeros(2*n,1);
    for ii=1:n
        if(Node(ii,2)==1) %slack
            y1(ii)=0;
            y2(ii)=0;
        end
        if(Node(ii,2)==2) %PV
            for jj=1:n
                y1(ii)=y1(ii)+e(ii)*(G(ii,jj)*e(jj)-B(ii,jj)*f(jj))+f(ii)*(G(ii,jj)*f(jj)+B(ii,jj)*e(jj));
                y2(ii)=e(ii)^2+f(ii)^2;
            end
        end
        if(Node(ii,2)==3) %PQ
           for jj=1:n
              y1(ii)= y1(ii)+e(ii)*(G(ii,jj)*e(jj)-B(ii,jj)*f(jj))+f(ii)*(G(ii,jj)*f(jj)+B(ii,jj)*e(jj));
              y2(ii)= y2(ii)+f(ii)*(G(ii,jj)*e(jj)-B(ii,jj)*f(jj))-e(ii)*(G(ii,jj)*f(jj)+B(ii,jj)*e(jj));
           end
        end
    end
    for ii=1:n
        y(2*ii-1)=y1(ii);
        y(2*ii)=y2(ii);
    end
    

    %build jacobian
    J=zeros(2*n);
    for ii=1:n
        if(Node(ii,2)==1) % slack
           for jj =1:n
               if(ii==jj)
                  J(2*ii-1,2*jj-1)=1;
                  J(2*ii,2*jj)=1;
               end
           end
        end
        if(Node(ii,2)==2)% PV
           for jj=1:n
              if(ii==jj)
                 for t=1:n
                    J(2*ii-1,2*jj-1)= J(2*ii-1,2*jj-1)+G(ii,t)*e(t)-B(ii,t)*f(t);
                    J(2*ii-1,2*jj)= J(2*ii-1,2*jj)+B(ii,t)*e(t)+G(ii,t)*f(t);
                 end
                 J(2*ii-1,2*jj-1)= J(2*ii-1,2*jj-1)+G(ii,jj)*e(ii)+B(ii,jj)*f(ii);
                 J(2*ii-1,2*jj)= J(2*ii-1,2*jj)-B(ii,jj)*e(ii)+G(ii,jj)*f(ii);
                 J(2*ii,2*jj-1)= 2*e(jj);
                 J(2*ii,2*jj)=2*f(jj); 
              end
              if(ii~=jj)
                  J(2*ii-1,2*jj-1)= G(ii,jj)*e(ii)+B(ii,jj)*f(ii);
                  J(2*ii-1,2*jj) = -B(ii,jj)*e(ii)+G(ii,jj)*f(ii);
                  J(2*ii,2*jj-1)= 0;
                  J(2*ii,2*jj)= 0; 
              end
           end
        end
        if(Node(ii,2)==3) %PQ
           for jj=1:n
              if(ii==jj)
                  for t=1:n
                     J(2*ii-1,2*jj-1)= J(2*ii-1,2*jj-1)+G(ii,t)*e(t)-B(ii,t)*f(t);
                     J(2*ii-1,2*jj)= J(2*ii-1,2*jj)+B(ii,t)*e(t)+G(ii,t)*f(t);
                     J(2*ii,2*jj-1)= J(2*ii,2*jj-1)-(B(ii,t)*e(t)+G(ii,t)*f(t));
                     J(2*ii,2*jj)= J(2*ii,2*jj)+G(ii,t)*e(t)-B(ii,t)*f(t);
                  end
                  J(2*ii-1,2*jj-1)= J(2*ii-1,2*jj-1)+G(ii,jj)*e(ii)+ B(ii,jj)*f(ii);
                  J(2*ii-1,2*jj)= J(2*ii-1,2*jj)-B(ii,jj)*e(ii)+ G(ii,jj)*f(ii);
                  J(2*ii,2*jj-1)= J(2*ii,2*jj-1)-B(ii,jj)*e(ii)+ G(ii,jj)*f(ii);
                  J(2*ii,2*jj)= J(2*ii,2*jj)-(G(ii,jj)*e(ii)+ B(ii,jj)*f(ii));
              end
              if(ii~=jj)
                  J(2*ii-1,2*jj-1)= G(ii,jj)*e(ii)+ B(ii,jj)*f(ii);
                  J(2*ii-1,2*jj)= -B(ii,jj)*e(ii)+ G(ii,jj)*f(ii);
                  J(2*ii,2*jj-1)= -B(ii,jj)*e(ii)+ G(ii,jj)*f(ii);
                  J(2*ii,2*jj)= -(G(ii,jj)*e(ii)+ B(ii,jj)*f(ii));
              end
          end
        end
    end
    disp(J);
    %to get a new x
    deltay=y-ys;
    disp(deltay);
    deltax=J^(-1)*deltay;
    disp(deltax);
    x=x-deltax;
    %new e f
    for a=1:n
        e(a)=x(2*a-1);
        f(a)=x(2*a);
    end

    %loop judgement
    max=0;
    precision=10^(-6);
    for a=1:2*n
        if(abs(deltax(a))>=max)
            max=abs(deltax(a));
        end
    end
    if (max<=precision)
        break;
    end
end
%result [number P Q V theta]
result=zeros(n,4);

for ii=1:n
     if(Node(ii,2)==1) %slack
         result(ii,1)=ii;
         for jj=1:n
            result(ii,2)= result(ii,2)+e(ii)*(G(ii,jj)*e(jj)-B(ii,jj)*f(jj))+f(ii)*(G(ii,jj)*f(jj)+B(ii,jj)*e(jj));
            result(ii,3)= result(ii,3)+f(ii)*(G(ii,jj)*e(jj)-B(ii,jj)*f(jj))-e(ii)*(G(ii,jj)*f(jj)+B(ii,jj)*e(jj));
         end
         result(ii,4)=Node(ii,5);
         result(ii,5)=atan(f(ii)/e(ii));
     end
     if(Node(ii,2)==2)% PV
         result(ii,1)=ii;
         result(ii,2)= Node(ii,3);
         for jj=1:n
            result(ii,3)= result(ii,3)+f(ii)*(G(ii,jj)*e(jj)-B(ii,jj)*f(jj))-e(ii)*(G(ii,jj)*f(jj)+B(ii,jj)*e(jj));
         end
         result(ii,4)= Node(ii,5);
         result(ii,5)= atan(f(ii)/e(ii));
     end
     if(Node(ii,2)==3)% PQ
         result(ii,1)=ii;
         result(ii,2)= Node(ii,3);
         result(ii,3)= Node(ii,4);
         result(ii,4)= (e(ii)^2+f(ii)^2)^0.5;
         result(ii,5)= atan(f(ii)/e(ii));
     end
end
