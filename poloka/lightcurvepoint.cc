#include <iostream>
#include <poloka/lightcurvepoint.h>

LightCurvePoint::LightCurvePoint() {
  flux=0;
  eflux=0;
  modifiedjulianday=0;
  mag=99;
  emag_minus=0;
  emag_plus=0;
}
