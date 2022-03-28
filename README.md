# satcompet
Adriana Sanz 22/05/2017  
Viral competition with a satellite using a simple compartmental model of two virus in competition, one of them assisted by a satellite.
<dl>
<dt>Fist order system of equations:</dt>
	<dd>h' = g - (d +p_x*x+p_y*y + p_sy*s)*h</dd>
	<dd>x' = (p_x*h - (d + d_x))*x</dd>
	<dd>y' = (p_y*h - (d + d_y) - p_s*s)*y + p_Y*s*h </dd>
	<dd>s' = (p_s*y -(d + d_s) + p_sy*h)*s</dd>
<dt>Variables:</dt>
	<dd>h heathy hosts</dd>
	<dd>x infected hosts by virus x (competitor virus)</dd> 
	<dd>y infected hosts by virus y (helper virus)</dd>
	<dd>s infected hosts by virus y and satellite</dd> 
<dt>Parameters:</dt>
	<dd>g income rate of uninfected hosts</dd> 
	<dd>p_i transmision rate i &#8712; {x,y,s}</dd> 
	<dd>d_i decay rate i &#8712; {x,y,s}</dd>   
</dl>

### 1. satcompet.c
Implements Runge Kutta integration step of the system. Parameters used are: d=0.05; p_x=0.1; p_y=0.1; d_y=0.25; p_s=0.8; p_{sy}=0.04; d_s=0.3. We study three phenomena:    
  - 1.A - Commensalism: the satellite does not affect the viral competition - use p_Y = 0.0766  
  - 1.B - Mutualism - bistability: the satellite helps the helper virus in the viral competition - use p_Y=0.2  
  - 1.C - Parasitism - coexistence: the satellite helps the competitor virus in the viral competition - use p_Y=0.16
	
### 2. Compilation routine: 
```
$ gcc -std=c99 satcompet.c -o nameofile
```
### 3. Run (WARNING:this could take ~ 1 day to run)
```
$ ./nameofile
```
### 4. Plot the results on matlab with fig-satellite.m
