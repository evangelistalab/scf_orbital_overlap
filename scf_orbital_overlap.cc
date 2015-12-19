/*
 *@BEGIN LICENSE
 *
 * scf_orbital_overlap by Psi4 Developer, a plugin to:
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>
#include <libscf_solver/rohf.h>

using namespace boost;

namespace psi{ namespace scf_orbital_overlap {

extern "C"
int read_options(std::string name, Options& options)
{
    if (name == "SCF_ORBITAL_OVERLAP"|| options.read_globals()) {
        /*- The smaller basis -*/
        options.add_str("BASIS1","");
        /*- The larger basis -*/
        options.add_str("BASIS2","");
        /*- The DF basis for basis 1 -*/
        options.add_str("DF_BASIS1", "");
        /*- The DF basis for basis 2 -*/
        options.add_str("DF_BASIS2", "");
    }

    return true;
}

extern "C"
PsiReturnType scf_orbital_overlap(Options& options)
{
    boost::shared_ptr<PSIO> psio = PSIO::shared_object();

    // Run SCF with basis #1
    options.set_global_str("BASIS",options.get_str("BASIS1"));
    options.set_global_str("DF_BASIS_SCF", options.get_str("DF_BASIS1"));
    outfile->Printf("\n  Switched to %s basis",options.get_str("BASIS").c_str());

    boost::shared_ptr<Wavefunction> scf_basis1(new scf::ROHF(options, psio));
    Process::environment.set_wavefunction(scf_basis1);
    double scf_basis1_energy = scf_basis1->compute_energy();
    boost::shared_ptr<Matrix> Ca1 = scf_basis1->Ca();

    // Run SCF with basis #2
    options.set_global_str("BASIS",options.get_str("BASIS2"));
    options.set_global_str("DF_BASIS_SCF", options.get_str("DF_BASIS2"));
    outfile->Printf("\n  Switched to %s basis",options.get_str("BASIS").c_str());

    PSIOManager::shared_object()->psiclean();

    boost::shared_ptr<Wavefunction> scf_basis2(new scf::ROHF(options, psio));
    Process::environment.wavefunction().reset();
    Process::environment.set_wavefunction(scf_basis2);
    double scf_basis2_energy = scf_basis2->compute_energy();
    boost::shared_ptr<Matrix> Ca2 = scf_basis2->Ca();

    // Compute the orbital overlap
    int nirrep = Process::environment.wavefunction()->nirrep();
    boost::shared_ptr<Molecule> molecule = Process::environment.molecule();
    boost::shared_ptr<BasisSet> Basis1 = BasisSet::pyconstruct_orbital(molecule, "BASIS",options.get_str("BASIS1"));
    boost::shared_ptr<BasisSet> Basis2 = BasisSet::pyconstruct_orbital(molecule, "BASIS",options.get_str("BASIS2"));

    boost::shared_ptr<IntegralFactory> integral_12(new IntegralFactory(Basis1, Basis2, Basis2, Basis2));
    boost::shared_ptr<IntegralFactory> integral_basis1(new IntegralFactory(Basis1, Basis1, Basis1, Basis1));
    boost::shared_ptr<IntegralFactory> integral_basis2(new IntegralFactory(Basis2, Basis2, Basis2, Basis2));

    boost::shared_ptr<SOBasisSet> soBasis1(new SOBasisSet(Basis1, integral_basis1));
    boost::shared_ptr<SOBasisSet> soBasis2(new SOBasisSet(Basis2, integral_basis2));

    Dimension nsopi_Basis1 = soBasis1->dimension();
    Dimension nsopi_Basis2 = soBasis2->dimension();

    boost::shared_ptr<MatrixFactory> soFactoryMixed(new MatrixFactory);
    soFactoryMixed->init_with(nsopi_Basis1,nsopi_Basis2);

    // Form the overlap matrix in the mixed basis
    boost::shared_ptr<OneBodySOInt> sOBI_cu(integral_12->so_overlap());
    SharedMatrix S12so(soFactoryMixed->create_matrix("Overlap"));
    sOBI_cu->compute(S12so);

    boost::shared_ptr<Matrix> overlap = Matrix::triplet(Ca1,S12so,Ca2,true,false,false);

    CharacterTable ct = Process::environment.molecule()->point_group()->char_table();


    outfile->Printf("\n  Orbital overlaps:");
    outfile->Printf("\n  --------------------------");
    outfile->Printf("\n   Basis 1    Basis 2 ");
    outfile->Printf("\n  --------------------------");
    for (int h = 0; h < nirrep; ++h){
        for (int i = 0; i < overlap->rowspi(h); ++i){
            double max = 0.0;
            int maxj = 0;
            for (int j = 0; j < overlap->colspi(h); ++j){
                if (std::fabs(overlap->get(h,i,j)) > max){
                    max = std::fabs(overlap->get(h,i,j));
                    maxj = j;
                }
            }
            outfile->Printf("\n  %3d-%-3s -> %3d-%-3s (%.2f)",i,
                            ct.gamma(h).symbol(),maxj,ct.gamma(h).symbol(),max * 100.0);
        }
    }
    outfile->Printf("\n  --------------------------");

    return Success;
}

}} // End namespaces

