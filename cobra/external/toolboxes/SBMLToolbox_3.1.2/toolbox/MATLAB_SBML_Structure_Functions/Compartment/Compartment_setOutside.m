function SBMLCompartment = Compartment_setOutside(SBMLCompartment, outside)
%
%   Compartment_setOutside 
%             takes  1) an SBMLCompartment structure 
%             and    2) a string representing the outside to be set
%
%             and returns 
%               the compartment with the outside set
%
%       SBMLCompartment = Compartment_setOutside(SBMLCompartment, 'outside')

%  Filename    :   Compartment_setOutside.m
%  Description :
%  Author(s)   :   SBML Development Group <sbml-team@caltech.edu>
%  $Id: Compartment_setOutside.m 7155 2008-06-26 20:24:00Z mhucka $
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
if (~isstruct(SBMLCompartment))
    error(sprintf('%s\n%s', ...
      'Compartment_setOutside(SBMLCompartment)', ...
      'argument must be an SBML compartment structure'));
end;
 
[sbmlLevel, sbmlVersion] = GetLevelVersion(SBMLCompartment);

if (~isSBML_Compartment(SBMLCompartment, sbmlLevel, sbmlVersion))
    error(sprintf('%s\n%s', 'Compartment_setOutside(SBMLCompartment, outside)', 'first argument must be an SBML compartment structure'));
elseif (~ischar(outside))
    error(sprintf('Compartment_setOutside(SBMLCompartment, outside)\n%s', 'second argument must be a string representing the outside of the compartment'));
end;

SBMLCompartment.outside = outside;