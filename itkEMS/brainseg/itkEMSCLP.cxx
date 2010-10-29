/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkEMSCLP.cxx,v $
  Language:  C++
  Date:      $Date: 2009/03/18 14:28:08 $
  Version:   $Revision: 1.5 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notices for more information.
  =========================================================================*/

extern int itkEMSCLPPrimary(int argc, char * argv[]);

//main function built in itkEMSCLPPrimary.cxx so that testing only builds templates once.
int main(int argc, char * argv[])
{
  return itkEMSCLPPrimary(argc,argv);
}

