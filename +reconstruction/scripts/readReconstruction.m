

[filename,pathname] = uigetfile('*.ASC') ;
cbuf = sprintf('%s%s', pathname, filename);

[tree, coords, contours, name, path] = neurolucida_tree(cbuf) ;

