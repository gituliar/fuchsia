.PRECIOUS: %_norm.m %_fuchs.m

FUCHSIA = sage -python ../fuchsia.py

# Examples from R.Lee
#
lee_1_y.m: lee_1.m
	$(FUCHSIA) changevar -m $@ $< "(1+y^2)/(1-y^2)"
lee_1_y_eps.m: lee_1_y.m
	$(FUCHSIA) reduce --use-maple -v -x y -m $@ -t $(basename $@)_t.m $<

lee_2_y.m: lee_2.m
	$(FUCHSIA) changevar -m $@ $< "(1-y)^2/(y)"
lee_2_y_eps.m: lee_2_y.m
	$(FUCHSIA) reduce --use-maple -v -x y -m $@ -t $(basename $@)_t.m $<

lee_3_eps.m: lee_3.m
	$(FUCHSIA) reduce --use-maple -m $@ -t $(basename $@)_t.m $<

# Example from C.Papadopolous
#
pap_1_eps.m: pap_1.m
	$(FUCHSIA) reduce --use-maple -e ep -m $@ -t $(basename $@)_t.m $<

# Other examples
#
%_eps.m: %.m
	$(FUCHSIA) reduce -m $@ -t $(basename $@)_t.m $<
%_fuchs.m: %.m
	$(FUCHSIA) fuchsify -m $@ -t $(basename $@)_t.m $<
%_norm.m: %_fuchs.m
	$(FUCHSIA) normalize -m $@ -t $(basename $@)_t.m $<
%_fact.m: %_norm.m
	$(FUCHSIA) factorize -m $@ -t $(basename $@)_t.m $<
