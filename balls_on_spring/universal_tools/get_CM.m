function [ x_cm, y_cm, z_cm ] = get_CM( m, x, varargin )
%Locates center of mass given row vectors of the particles x-, y- and z-components

L=size(x,1);
m_matrix=repmat(m, L,1);
M=sum(m);

switch nargin
    case 2
        x_cm=sum(m_matrix.*x, 2)/M;
        y_cm=[];z_cm=[];
    case 3
        y=varargin{1};
        x_cm=sum(m_matrix.*x, 2)/M;
        y_cm=sum(m_matrix.*y, 2)/M;
        z_cm=[];
    case 4
        y=varargin{1}; z=varargin{2};
        x_cm=sum(m_matrix.*x, 2)/M;
        y_cm=sum(m_matrix.*y, 2)/M;
        z_cm=sum(m_matrix.*z, 2)/M;
    otherwise
        x_cm=[];y_cm=[];z_cm=[];
        fprintf('Error: %d input arguments, max 4. \n\n', nargin);
end


end

