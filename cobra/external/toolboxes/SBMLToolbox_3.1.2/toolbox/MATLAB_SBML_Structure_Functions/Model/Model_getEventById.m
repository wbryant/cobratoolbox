function event = Model_getEventById(SBMLModel, id)
%
%   Model_getEventById 
%             takes  1) an SBMLModel structure 
%             and    2) a string representing the id of the event to be found
%
%             and returns 
%               the event structure with the matching id 
%               or an empty structure if no such event exists
%               
%       event = Model_getEventById(SBMLModel, 'id')

%  Filename    :   Model_getEventById.m
%  Description :
%  Author(s)   :   SBML Development Group <sbml-team@caltech.edu>
%  $Id: Model_getEventById.m 7155 2008-06-26 20:24:00Z mhucka $
%  $Source v $
%
%<!---------------------------------------------------------------------------
% This file is part of SBMLToolbox.  Please visit http://sbml.org for more
% information about SBML, and the latest version of SBMLToolbox.
%
% Copyright 2005-2007 California Institute of Technology.
% Copyright 2002-2005 California Institute of Technology and
%                     Japan Science and Technology Corporation.
% 
% This library is free software; you can redistribute it and/or modify it
% under the terms of the GNU Lesser General Public License as published by
% the Free Software Foundation.  A copy of the license agreement is provided
% in the file named "LICENSE.txt" included with this software distribution.
% and also available online as http://sbml.org/software/sbmltoolbox/license.html
%----------------------------------------------------------------------- -->



% check that input is correct
if (~isSBML_Model(SBMLModel))
    error(sprintf('%s\n%s', 'Model_getEventById(SBMLModel, id)', 'first argument must be an SBML model structure'));
elseif (~ischar(id))
    error(sprintf('%s\n%s', 'Model_getEventById(SBMLModel, id)', 'second argument must be a string'));
elseif (SBMLModel.SBML_level ~= 2)
    error(sprintf('%s\n%s', 'Model_getEventById(SBMLModel, id)', 'no events in a level 1 model'));   
end;

event = [];

for i = 1:length(SBMLModel.event)
    if (strcmp(id, SBMLModel.event(i).id))
        event = SBMLModel.event(i);
        return;
    end;
end;