function varargout = nrbderiv (nurbs)
% 
% NRBDERIV: Construct the first and second derivative representation of a
%           NURBS curve, surface or volume.
% 
% Calling Sequence:
% 
%   ders = nrbderiv (nrb);
%   [ders, ders2] = nrbderiv (nrb);
%   [ders, ders2, ders3] = nrbderiv (nrb);
% 
% INPUT:
% 
%   nrb		: NURBS data structure, see nrbmak.
%
% OUTPUT:
% 
%   ders:  A data structure that represents the first
% 		    derivatives of a NURBS curve, surface or volume.
%   ders2: A data structure that represents the second
% 		    derivatives of a NURBS curve, surface or volume.
%   ders3: A data structure that represents the third
% 		    derivatives of a NURBS curve, surface.
% 
% Description:
% 
%   The derivatives of a B-Spline are themselves a B-Spline of lower degree,
%   giving an efficient means of evaluating multiple derivatives. However,
%   although the same approach can be applied to NURBS, the situation for
%   NURBS is more complex. We have followed in this function the same idea
%   that was already used for the first derivative in the function nrbderiv.
%   The second derivative data structure can be evaluated later with the
%   function nrbdeval.
% 
% See also:
% 
%       nrbdeval
%
% Copyright (C) 2000 Mark Spink
% Copyright (C) 2010 Carlo de Falco
% Copyright (C) 2010, 2011 Rafael Vazquez
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 2 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

if (~isstruct(nurbs))
  error('NURBS representation is not structure!');
end

if (~strcmp(nurbs.form,'B-NURBS'))
  error('Not a recognised NURBS representation');
end

% We raise the degree to avoid errors in the computation of the second
% derivative
if (iscell (nurbs.knots))
  ndim = size(nurbs.knots, 2);
else
  ndim = 1;
end

if (nargout == 2)
  degelev  = max (2*ones(1, ndim) - (nurbs.order-1), 0);
  nurbs    = nrbdegelev (nurbs, degelev);
elseif (nargout == 3)
  degelev  = max (3*ones(1, ndim) - (nurbs.order-1), 0);
  nurbs    = nrbdegelev (nurbs, degelev);
end

degree = nurbs.order - 1;

keyboard


if (ndim == 3)
%% NURBS structure represents a volume
  num1 = nurbs.number(1);
  num2 = nurbs.number(2);
  num3 = nurbs.number(3);

% taking derivatives along the u direction
  dknots = nurbs.knots;
  dcoefs = permute (nurbs.coefs,[1 3 4 2]);
  dcoefs = reshape (dcoefs,4*num2*num3,num1);
  [dcoefs,dknots{1}] = bspderiv (degree(1),dcoefs,nurbs.knots{1});
  dcoefs = permute (reshape (dcoefs,[4 num2 num3 size(dcoefs,2)]),[1 4 2 3]);
  dnurbs{1} = nrbmak (dcoefs, dknots);

  if (nargout == 2)
% taking second derivative along the u direction (duu)
    dknots2 = dknots;
    dcoefs2 = permute (dcoefs, [1 3 4 2]);
    dcoefs2 = reshape (dcoefs2, 4*num2*num3, []);
    [dcoefs2, dknots2{1}] = bspderiv (degree(1)-1, dcoefs2, dknots{1});
    dcoefs2 = permute (reshape (dcoefs2, 4, num2, num3, []), [1 4 2 3]);
    dnurbs2{1,1} = nrbmak (dcoefs2, dknots2); 

% taking second derivative along the v direction (duv and dvu)
    dknots2 = dknots;
    dcoefs2 = permute (dcoefs,[1 2 4 3]);
    dcoefs2 = reshape (dcoefs2, 4*(num1-1)*num3, num2);
    [dcoefs2, dknots2{2}] = bspderiv (degree(2), dcoefs2, dknots{2});
    dcoefs2 = permute (reshape (dcoefs2, 4, num1-1, num3, []), [1 2 4 3]);
    dnurbs2{1,2} = nrbmak (dcoefs2, dknots2);
    dnurbs2{2,1} = dnurbs2{1,2};

% taking second derivative along the w direction (duw and dwu)
    dknots2 = dknots;
    dcoefs2 = reshape (dcoefs, 4*(num1-1)*num2, num3);
    [dcoefs2, dknots2{3}] = bspderiv (degree(3), dcoefs2, dknots{3});
    dcoefs2 = reshape (dcoefs2, 4, num1-1, num2, []);
    dnurbs2{1,3} = nrbmak (dcoefs2, dknots2);
    dnurbs2{3,1} = dnurbs2{1,3};
  end

% taking derivatives along the v direction
  dknots = nurbs.knots;
  dcoefs = permute (nurbs.coefs,[1 2 4 3]);
  dcoefs = reshape (dcoefs,4*num1*num3,num2);
  [dcoefs,dknots{2}] = bspderiv (degree(2),dcoefs,nurbs.knots{2});
  dcoefs = permute (reshape (dcoefs,[4 num1 num3 size(dcoefs,2)]),[1 2 4 3]);
  dnurbs{2} = nrbmak (dcoefs, dknots);

  if (nargout == 2)
% taking second derivative along the v direction (dvv)
    dknots2 = dknots;
    dcoefs2 = permute (dcoefs,[1 2 4 3]);
    dcoefs2 = reshape (dcoefs2, 4*num1*num3, num2-1);
    [dcoefs2, dknots2{2}] = bspderiv (degree(2)-1, dcoefs2, dknots{2});
    dcoefs2 = permute (reshape (dcoefs2, 4, num1, num3, []), [1 2 4 3]);
    dnurbs2{2,2} = nrbmak (dcoefs2, dknots2);

% taking second derivative along the w direction (dvw and dwv)
    dknots2 = dknots;
    dcoefs2 = reshape (dcoefs, 4*num1*(num2-1), num3);
    [dcoefs2, dknots2{3}] = bspderiv (degree(3), dcoefs2, dknots{3});
    dcoefs2 = reshape (dcoefs2, 4, num1, num2-1, []);
    dnurbs2{2,3} = nrbmak (dcoefs2, dknots2);
    dnurbs2{3,2} = dnurbs2{2,3};
  end

% taking derivatives along the w direction
  dknots = nurbs.knots;
  dcoefs = reshape (nurbs.coefs,4*num1*num2,num3);
  [dcoefs,dknots{3}] = bspderiv (degree(3),dcoefs,nurbs.knots{3});
  dcoefs = reshape (dcoefs,[4 num1 num2 size(dcoefs,2)]);
  dnurbs{3} = nrbmak (dcoefs, dknots);

  if (nargout == 2)
% taking second derivative along the w direction (dww)
    dknots2 = dknots;
    dcoefs2 = reshape (dcoefs, 4*num1*num2, num3-1);
    [dcoefs2, dknots2{3}] = bspderiv (degree(3)-1, dcoefs2, dknots{3});
    dcoefs2 = reshape (dcoefs2, 4, num1, num2, []);
    dnurbs2{3,3} = nrbmak (dcoefs2, dknots2);
  end

elseif (ndim == 2)

  %% NURBS structure represents a surface
  num1 = nurbs.number(1);
  num2 = nurbs.number(2);

%% R_uuu
% taking first derivative along the u direction
  dknots = nurbs.knots;
  dcoefs = permute (nurbs.coefs,[1 3 2]);
  dcoefs = reshape (dcoefs,4*num2,num1);
  [dcoefs,dknots{1}] = bspderiv (degree(1),dcoefs,nurbs.knots{1});
  dcoefs = permute (reshape (dcoefs,[4 num2 size(dcoefs,2)]),[1 3 2]);
  dnurbs{1} = nrbmak (dcoefs, dknots);

  if (nargout >= 2)
% taking second derivative along the u direction (duu)
    dknots2 = dknots;
    
    dcoefs2 = permute (dcoefs, [1 3 2]);
    
    dcoefs2 = reshape (dcoefs2, 4*num2, num1-1);
    
    [dcoefs2, dknots2{1}] = bspderiv (degree(1)-1, dcoefs2, dknots{1});
    
    dcoefs2 = permute (reshape (dcoefs2, 4, num2, []), [1 3 2]);
    
    dnurbs2{1,1} = nrbmak (dcoefs2, dknots2); 
        
  if (nargout >= 3)
      % taking third derivative along the u direction (duuv)
    dknots3 = dknots2;
    dcoefs3 = permute (dcoefs3, [1 3 2]);
    dcoefs3 = reshape (dcoefs3, 4*(num2-1), num1-1);
    [dcoefs3, dknots3{1}] = bspderiv (degree(1)-1, dcoefs3, dknots2{2});
    dcoefs3 = reshape (dcoefs3, 4, num1-2, []);
    dnurbs3{1,2} = nrbmak (dcoefs3, dknots3);
   
% taking third derivative along the u direction (duuu) 
    dknots3 = dknots2;
    dcoefs3 = permute (dcoefs2, [1 3 2]);
    dcoefs3 = reshape (dcoefs3, 4*num2, []);
    [dcoefs3, dknots3{1}] = bspderiv (degree(1)-2, dcoefs3, dknots2{1});
    dcoefs3 = permute (reshape (dcoefs3, 4, num2, []), [1 3 2]);
    dnurbs3{1,1} = nrbmak (dcoefs3, dknots3);
  end

% taking second derivative along the v direction (duv and dvu)
    dknots2 = dknots;
    dcoefs2 = reshape (dcoefs, 4*(num1-1), num2);
    [dcoefs2, dknots2{2}] = bspderiv (degree(2), dcoefs2, dknots{2});
    dcoefs2 = reshape (dcoefs2, 4, num1-1, []);
    dnurbs2{1,2} = nrbmak (dcoefs2, dknots2);
    dnurbs2{2,1} = dnurbs2{1,2};
  end

% taking first derivative along the v direction
  dknots = nurbs.knots;
  dcoefs = reshape (nurbs.coefs,4*num1,num2);
  [dcoefs,dknots{2}] = bspderiv (degree(2),dcoefs,nurbs.knots{2});
  dcoefs = reshape (dcoefs,[4 num1 size(dcoefs,2)]);
  dnurbs{2} = nrbmak (dcoefs, dknots);

  if (nargout >= 2)
% taking second derivative along the v direction (dvv)
    dknots2 = dknots;
    dcoefs2 = reshape (dcoefs, 4*num1, num2-1);
    [dcoefs2, dknots2{2}] = bspderiv (degree(2)-1, dcoefs2, dknots{2});
    dcoefs2 = reshape (dcoefs2, 4, num1, []);
    dnurbs2{2,2} = nrbmak (dcoefs2, dknots2);
  end
  
  if (nargout >= 3)
      % taking third derivative along the v direction (duvv)    
    dknots3 = dknots2;
    dcoefs3 = reshape (dcoefs2, 4*(num1-1), num2-1);
    [dcoefs3, dknots3{2}] = bspderiv (degree(2)-1, dcoefs3, dknots2{2});
    dcoefs3 = reshape (dcoefs3, 4, num1-1, []);
    dnurbs3{2,1} = nrbmak (dcoefs3, dknots3);
% taking third derivative along the v direction (dvvv)
      dknots3 = dknots2;
      dcoefs3 = reshape (dcoefs2, 4*num1, num2-2);
      [dcoefs3, dknots3{2}] = bspderiv (degree(2)-2, dcoefs3, dknots2{2});
      dcoefs3 = reshape (dcoefs3, 4, num1, []);
      dnurbs3{2,2} = nrbmak (dcoefs3, dknots3);
  end
%%  
else
%% NURBS structure represents a curve
  [dcoefs,dknots] = bspderiv (degree, nurbs.coefs, nurbs.knots);
  dnurbs = nrbmak (dcoefs, dknots);
  if (nargout >= 2)
    [dcoefs2,dknots2] = bspderiv (degree-1, dcoefs, dknots);
    dnurbs2 = nrbmak (dcoefs2, dknots2);
  end
  if (nargout == 3)
    [dcoefs3,dknots3] = bspderiv (degree-2, dcoefs2, dknots2);
    dnurbs3 = nrbmak (dcoefs3, dknots3);
  end
end

varargout{1} = dnurbs;
if (nargout == 2)
  varargout{2} = dnurbs2;
end
if (nargout == 3)
    varargout{2} = dnurbs2;
    varargout{3} = dnurbs3;
end

end