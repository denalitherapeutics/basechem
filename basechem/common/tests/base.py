import glob
import json
import os
import shutil
import tempfile
import unittest.mock as mock
from io import BytesIO

from django.conf import settings
from django.core.files.uploadedfile import SimpleUploadedFile
from django.test import TestCase, override_settings
from django.urls import reverse
from rdkit import Chem

from basechem.common.mmpdb_utils import initialize_mmpdb
from basechem.common.mocks.mock_inductive_utils import (
    mock_run_inductive_alogd_predict,
    mock_run_inductive_lm_predict,
    mock_update_inductive_logd_data,
)
from basechem.common.mocks.mock_mayachem_utils import mock_generate_torsion_alerts_NOOP
from basechem.common.tests.test_constants import BAD_MOLTEXT, MMP_MOLTEXT, MOLTEXT
from basechem.main.constants import DTX_PROPERTIES_TO_SAVE, MMP
from basechem.main.models.project_models import Project
from basechem.main.tests.factories import (
    CollectionFactory,
    CompoundFactory,
    ProjectFactory,
    SeriesFactory,
)
from basechem.users.tests.factories import BasechemUserFactory


class BasechemViewTestMixin:
    def setUp(self):
        super().setUp()
        self.user = BasechemUserFactory(
            username="test@example.com", email="test@example.com"
        )
        self.user.set_password("testpassword")
        self.user.save()
        self.client.login(username=self.user.username, password="testpassword")
        Project.objects.create(code="TEST")
        # for an unknown reason switching this to .get is not working
        self.test_project = Project.objects.filter(code="TEST")[0].pk

    def upload_sdf(self, url, metadata=None):
        """
        Creates a file and posts it to provided url view
        :param url: url to upload to
        :param metadata: additional data to add to the form upload
        """
        file_contents = open("basechem/main/tests/testdata/test_onecomp.sdf", "rb")
        test_sdf_file = SimpleUploadedFile("test_file.sdf", file_contents.read())
        file_data = {"upload_file": test_sdf_file}

        data = {
            "project": self.test_project,
            "upload_file": test_sdf_file,
        }

        if metadata:
            for k, v in metadata.items():
                data[k] = v

        response = self.client.post(url, data=data, files=file_data, follow=True)
        return response

    def upload_bad_sdf(self, url):
        """
        Creates a file and posts it to provided url view
        """
        file_contents = BytesIO(b"some dummy bcode data: \x00\x01")
        test_sdf_file = SimpleUploadedFile("test_file.sdf", file_contents.read())
        file_data = {"upload_file": test_sdf_file}

        data = {
            "project": self.test_project,
            "upload_file": test_sdf_file,
        }

        response = self.client.post(url, data=data, files=file_data, follow=True)
        return response

    def upload_mol(self, url, data=None):
        """
        Creates a mol text and posts it to submit compounds view
        :param url: url to upload to
        :param data: additional data to add to the form upload
        """
        form_data = {"project": self.test_project, "moltext": MOLTEXT}
        if data:
            form_data.update(data)
        response = self.client.post(url, data=form_data, follow=True)
        return response

    def upload_mmp_mol(self, variable={}):
        """
        Posts a structure as moltext and a variable region as a list of selected atom indices
        to the MMP SubmitCompoundsView
        :param variable: dict with the key "atoms" corresponding to a list of the selected variable region atom indices for the MMP search
        :return: the response of post to MMP submit compounds view
        """
        url = reverse("submit", kwargs={"nextview": MMP})
        index_map = {i: i for i in range(8)}
        sketcher_data = {
            "moltext": MMP_MOLTEXT,
            "variable": variable,
            "maps": {"atoms": index_map, "bonds": index_map},
        }
        data = {"project": self.test_project, "sketcher": json.dumps(sketcher_data)}
        response = self.client.post(url, data=data, follow=True)
        return response

    def upload_bad_mol(self, url):
        """
        Creates a mol text and posts it to SubmitCompoundsView
        """
        data = {"project": self.test_project, "moltext": BAD_MOLTEXT}
        response = self.client.post(url, data=data, follow=True)
        return response


class BasechemFormTestMixin:
    def setUp(self):
        super().setUp()
        Project.objects.create(code="TEST")
        self.test_project = Project.objects.all()[0].pk

        # GOOD SDF #
        file_contents = open("basechem/main/tests/testdata/test_onecomp.sdf", "rb")
        self.test_sdf_file = SimpleUploadedFile("test_file.sdf", file_contents.read())
        self.file_data = {"upload_file": self.test_sdf_file}

        # BAD SDF #
        file_contents = BytesIO(b"some dummy bcode data: \x00\x01")
        self.test_bad_sdf = SimpleUploadedFile(
            "test_bad_file.sdf", file_contents.read()
        )
        self.bad_sdf_data = {"upload_file": self.test_bad_sdf}

        # BAD OTHER FILE #
        file_contents = BytesIO(b"some dummy bcode data: \x00\x01")
        self.test_bad_file = SimpleUploadedFile(
            "test_bad_file.csv", file_contents.read()
        )
        self.bad_file_data = {"upload_file": self.test_bad_file}

        # MOL #
        self.bad_moltext = BAD_MOLTEXT
        self.moltext = MOLTEXT

        self.user = BasechemUserFactory(
            username="test@example.com", email="test@example.com"
        )


@override_settings(MEDIA_ROOT=tempfile.mkdtemp(prefix="{}/tmp".format(os.getcwd())))
@override_settings(
    SLURM_SHARED_FILES_TMP_DIR="/efs-mount/groups/compchem/shared_files/basechem_tmp/local_testsuite"
)
class BasechemNoMockTestCase(TestCase):
    """
    This test case should only be used for testing utils (ex InductiveBio), as it
    does not include any mocked functions. Most tests (and all tests in the `main` app) should use
    the `BasechemTestCase` instead of this.
    """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()

        cls.project1 = ProjectFactory(code="PyTest")
        cls.admin_user = BasechemUserFactory(
            first_name="ADMIN",
            last_name="ADMIN",
            username=settings.ADMIN_USER,
            password=settings.ADMIN_PASSWORD,
            email=settings.ADMIN_EMAIL,
        )

        cls.owner = BasechemUserFactory()
        cls.collection = CollectionFactory(owner=cls.owner, project=cls.project1)
        cls.collection_no_cpds = CollectionFactory(
            owner=cls.owner, project=cls.project1
        )
        cls.collection_one_cpd = CollectionFactory(
            owner=cls.owner, project=cls.project1
        )
        cls.pub_collection = CollectionFactory(owner=cls.owner, project=cls.project1)
        cls.collection_2d = CollectionFactory(owner=cls.owner, project=cls.project1)
        cls.collection_3d = CollectionFactory(owner=cls.owner, project=cls.project1)
        cls.collection_torsion = CollectionFactory(
            owner=cls.owner, project=cls.project1
        )
        cls.collection_mmp = CollectionFactory(owner=cls.owner, project=cls.project1)

        one_comp_3d_file = open(
            "basechem/main/tests/testdata/test_onecomp_3d.sdf", "rb"
        )
        one_comp_2d_file = open(
            "basechem/main/tests/testdata/test_onecomp_2d.sdf", "rb"
        )
        one_comp_file = open("basechem/main/tests/testdata/test_onecomp.sdf", "rb")
        series_comp_file = open("basechem/main/tests/testdata/7eri_ligand.sdf", "rb")
        torsion_file = open("basechem/main/tests/testdata/test_torsion_input.sdf", "rb")
        three_comp_file = open(
            "basechem/main/tests/testdata/test_threecomp_mol.sdf", "rb"
        )
        receptor_file = open("basechem/main/tests/testdata/7eri_protein.pdb", "rb")
        receptor_file_mol2 = open(
            "basechem/main/tests/testdata/7eri_protein.mol2", "rb"
        )

        cls.one_3d_uploaded_file = SimpleUploadedFile(
            "one_comp.sdf", one_comp_3d_file.read()
        )
        cls.one_2d_uploaded_file = SimpleUploadedFile(
            "one_comp.sdf", one_comp_2d_file.read()
        )
        cls.one_uploaded_file = SimpleUploadedFile("one_comp.sdf", one_comp_file.read())
        cls.torsion_uploaded_file = SimpleUploadedFile(
            "one_comp.sdf", torsion_file.read()
        )
        cls.three_uploaded_file = SimpleUploadedFile(
            "three_comp.sdf", three_comp_file.read()
        )
        cls.receptor_uploaded_file = SimpleUploadedFile(
            "7eri_protein.pdb", receptor_file.read()
        )
        cls.receptor_mol2_uploaded_file = SimpleUploadedFile(
            "7eri_protein.mol2", receptor_file_mol2.read()
        )
        cls.series1_bound_file = SimpleUploadedFile(
            "7eri_ligand.sdf", series_comp_file.read()
        )

        # Prepare 3 mol Collection
        three_mols = cls.collection.handle_sdf_upload(cls.three_uploaded_file)
        cls.collection.handle_romols(three_mols, test=True)
        cls.collection.save()
        cls.compounds = list(cls.collection.compounds())

        # Prepare one mol Collection
        one_mol = cls.collection_one_cpd.handle_sdf_upload(cls.one_uploaded_file)
        cls.collection_one_cpd.handle_romols(one_mol, test=True)
        cls.collection_one_cpd.save()

        # Prepare 2D Collection
        one_mol_2d = cls.collection_2d.handle_sdf_upload(cls.one_2d_uploaded_file)
        cls.collection_2d.handle_romols(one_mol_2d, test=True)
        cls.collection_2d.save()

        # Prepare 3D Collection
        one_mol_3d = cls.collection_3d.handle_sdf_upload(cls.one_3d_uploaded_file)
        cls.collection_3d.handle_romols(one_mol_3d, test=True)
        cls.collection_3d.save()

        # Prepare Torsion Collection
        torsion_mol = cls.collection_torsion.handle_sdf_upload(
            cls.torsion_uploaded_file
        )
        cls.collection_torsion.handle_romols(torsion_mol, test=True)
        cls.collection_torsion.save()

        # Prepare MMP Collection
        mmp_mols = cls.create_assayed_compounds(cls)
        cls.collection_mmp.handle_romols(mmp_mols, test=True)
        mmp_comp = cls.collection_mmp.compounds().get(smiles="CCC1CCC(CF)CC1")
        cls.collection_mmp.metadata = {
            "mmp_analysis": {
                str(mmp_comp.pk): {
                    "constant_smiles": "*C1CCC(CC)CC1",
                    "variable_smiles": "*CF",
                }
            }
        }
        cls.collection_mmp.save()

        cls.series1 = SeriesFactory(
            dn_id="DN0000001",
            name="7ERI",
            smiles="ClC1=C([H])C(N=C(SC([H])([H])[H])N2[H])=C2C([H])=C1OC3=C(Cl)C(Cl)=C([H])C([H])=C3[H]",
            project=cls.project1,
            bound_state_file=cls.series1_bound_file,
            receptor_file=cls.receptor_uploaded_file,
            receptor_file_mol2=cls.receptor_mol2_uploaded_file,
            receptor_residues="15, 17, 108, 109, 110, 117, 118, 119, 201",
        )

        cls.series2 = SeriesFactory(
            dn_id="DN0000002",
            smiles="C1CCCCC1",
            project=cls.project1,
            ground_state_file=cls.one_uploaded_file,
        )

    def setUp(self):
        super().setUp()
        # Technically all the compounds *should have the same project but we don't specify one
        # so all these have None project for testing until we set them
        self.compound_a = self.compounds[0]
        self.compound_a.project = self.project1

        self.compound_b = self.compounds[1]
        self.compound_c = self.compounds[2]

    @classmethod
    def setUpClass(cls):
        """
        Set up mock for PUT requests, since we never want to use that in tests
        """
        super().setUpClass()
        initialize_mmpdb("test")
        # Set up mock patch for InductiveBio
        cls.main_update_inductive_logd_data_patch = mock.patch(
            "basechem.main.tasks.update_inductive_logd_data",
            mock_update_inductive_logd_data,
        )
        cls.main_update_inductive_logd_data_patch.start()

    @classmethod
    def tearDownClass(cls):
        super().tearDownClass()
        # Tear down mock patch for InductiveBio PUT
        cls.main_update_inductive_logd_data_patch.stop()
        for dirpath in glob.glob(f"{settings.PROJECT_ROOT}/tmp*"):
            shutil.rmtree(dirpath)

    def create_assayed_compounds(cls):
        """
        Creates new compounds from a test dataset with MMPs and returns the list of mols that should be used to
        create an MMP collection. Each compound in the collection has MMP compounds found in this test dataset.
        :returns: a list of mol objects that should be in the assay collection
        """
        test_sdf = "basechem/main/tests/testdata/test_mmps_w_data.sdf"
        assayed_mols = []
        suppl = Chem.SDMolSupplier(test_sdf)
        for mol in suppl:
            new_compound = CompoundFactory(
                inchi=Chem.MolToInchi(mol, options="/suu"),
                smiles=Chem.MolToSmiles(mol),
                dn_id=mol.GetProp("ID"),
                project=cls.project1,
            )
            for prop in DTX_PROPERTIES_TO_SAVE:
                if mol.HasProp(prop):
                    keyprop = prop.lower().replace(" ", "_")
                    new_compound.measured_data[keyprop] = mol.GetProp(prop)
            new_compound.save()
            if mol.GetProp("this week assayed") == "True":
                assayed_mols.append((mol, True))

        return assayed_mols


@override_settings(MEDIA_ROOT=tempfile.mkdtemp(prefix="{}/tmp".format(os.getcwd())))
@override_settings(SLURM_CONDA_EXEC_PATH="/usr/local/bin/conda")
@override_settings(
    SHARED_FILES_TMP_DIR_FROM_CONTAINER="/opt/shared_files/local_testsuite"
)
@override_settings(SLURM_REST_API_HOST="000")
class BasechemTestCase(BasechemNoMockTestCase):
    """
    Test case that should be used for all tests except those that directly test external
    code (InductiveBio, Mayachem, etc.). This expands on the `BasechemNoMockTestCase`
    by including patches functions that take a long time or can't be run on Github
    """

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        # Set up mock patch for InductiveBio
        cls.main_run_inductive_lm_predict_patch = mock.patch(
            "basechem.main.models.compound_models.run_inductive_lm_predict",
            mock_run_inductive_lm_predict,
        )
        cls.common_run_inductive_lm_predict_patch = mock.patch(
            "basechem.common.propcalc_utils.run_inductive_lm_predict",
            mock_run_inductive_lm_predict,
        )
        cls.main_run_inductive_alogd_predict_patch = mock.patch(
            "basechem.main.models.compound_models.run_inductive_alogd_predict",
            mock_run_inductive_alogd_predict,
        )
        cls.main_run_inductive_lm_predict_patch.start()
        cls.common_run_inductive_lm_predict_patch.start()
        cls.main_run_inductive_alogd_predict_patch.start()

        # Set up mock patch for Mayachem
        cls.generate_torsion_alerts_patch = mock.patch(
            "basechem.main.models.compound_models.generate_torsion_alerts",
            mock_generate_torsion_alerts_NOOP,
        )

        cls.generate_torsion_alerts_patch.start()

        # Set up mock patches for Dotmatics
        cls.collection_get_agg_ic50_data_patch = mock.patch(
            "basechem.main.models.collection_models.get_agg_ic50_data",
        )
        cls.collection_get_agg_ic50_data_mock = (
            cls.collection_get_agg_ic50_data_patch.start()
        )

    @classmethod
    def tearDownClass(cls):
        super().tearDownClass()
        # Tear down mock patch for InductiveBio
        cls.main_run_inductive_lm_predict_patch.stop()
        cls.common_run_inductive_lm_predict_patch.stop()
        cls.main_run_inductive_alogd_predict_patch.stop()

        # Tear down mock patch for Mayachem
        cls.generate_torsion_alerts_patch.stop()

        # Tear down mock patch for Dotmatics
        cls.collection_get_agg_ic50_data_patch.stop()
