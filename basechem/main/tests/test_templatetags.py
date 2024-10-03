import datetime

from django.core.files.uploadedfile import SimpleUploadedFile
from django.db.models import Prefetch
from django.template import Context, Template
from django.test import TestCase

from basechem.main.models.compound_models import Compound, CompoundOccurrence
from basechem.main.tests.factories import (
    BasechemUserFactory,
    CollectionFactory,
    CompoundFactory,
    ProjectFactory,
)


class TemplatetagTestCase(TestCase):
    def setUp(self):
        super().setUp()
        self.project = ProjectFactory()
        self.owner1 = BasechemUserFactory()
        self.owner2 = BasechemUserFactory()

        self.collection1 = CollectionFactory(owner=self.owner1, project=self.project)
        self.collection2 = CollectionFactory(owner=self.owner2, project=self.project)

        one_comp_file = open("basechem/main/tests/testdata/test_onecomp.sdf", "rb")
        one_uploaded_file = SimpleUploadedFile("one_comp.sdf", one_comp_file.read())
        one_mol = self.collection1.handle_sdf_upload(one_uploaded_file)

        self.collection1.handle_romols(one_mol, test=True)
        self.collection2.handle_romols(one_mol, test=True)
        self.collection1.save()
        self.collection2.save()

        # Sort compoundoccurrence_set by pk, because the `ProjectView`
        # does this and  `get_collection_urls`, `get_most_recent_date`, and `get_user_ids`
        # are only used in that view
        co_w_owner = CompoundOccurrence.objects.select_related("owner").order_by("pk")
        self.compounds = Compound.objects.all().prefetch_related(
            Prefetch("compoundoccurrence_set", queryset=co_w_owner)
        )

    def test_as_svg(self):
        """
        Test that the `as_svg` tag returns a compound as an svg
        """
        c = CompoundFactory(smiles="CC")
        out = Template("{% load tags %}" '{{ compound|as_svg:"200,50" }}').render(
            Context({"compound": c})
        )
        self.assertNotEqual(out, c.get_inline_svg(50, 50))  # wrong dimensions
        self.assertEqual(out, c.get_inline_svg(200, 50))  # right dimensions

    def test_as_transparent_svg(self):
        """
        Test that the `as_transparent_svg` tag returns a compound as an svg with a transparent background
        """
        c = CompoundFactory(smiles="CC")
        out = Template(
            "{% load tags %}" '{{ compound|as_transparent_svg:"200,50" }}'
        ).render(Context({"compound": c}))
        # wrong dimensions
        self.assertNotEqual(out, c.get_inline_svg(50, 50, transparent=True))
        # right dimensions
        self.assertEqual(out, c.get_inline_svg(200, 50, transparent=True))

    def test_get_collection_urls(self):
        """
        Test that the `get_collection_urls` tag returns all urls for collections which
        contain the compound
        """
        compound = self.compounds[0]
        out = Template("{% load tags %}" "{{ compound|get_collection_urls }}").render(
            Context({"compound": compound})
        )
        expected_out = f"<a href=/align/{self.collection1.id}>{self.collection1.id}</a>, <a href=/align/{self.collection2.id}>{self.collection2.id}</a>"
        self.assertEqual(out, expected_out)

    def test_get_most_recent_date(self):
        """
        Test that the `get_most_recent_date` tag returns most recent date
        of the compounds COs
        """
        compound = self.compounds[0]
        co = CompoundOccurrence.objects.filter(compound=compound)[0]
        co.generated = co.generated - datetime.timedelta(3)
        co.save()
        out = Template("{% load tags %}" "{{ compound|get_most_recent_date }}").render(
            Context({"compound": compound})
        )
        maxdate = max(
            [c.generated for c in CompoundOccurrence.objects.filter(compound=compound)]
        )
        expected_out = maxdate.strftime("%b. %d, %Y")
        self.assertEqual(out, expected_out)

    def test_get_user_ids(self):
        """
        Test that the `get_user_ids` tag returns all users with a COs associated
        with the compound
        """
        compound = self.compounds[0]
        cos = CompoundOccurrence.parents.filter(compound=compound)
        out = Template("{% load tags %}" "{{ compound|get_user_ids }}").render(
            Context({"compound": compound})
        )
        expected_out = ", ".join([co.owner.last_name for co in cos])
        self.assertEqual(sorted(out), sorted(expected_out))

    def test_get_delta_energy(self):
        """
        Test that the `get_delta_energy` tag returns the value of the delta_energy
        """
        data = {
            "co-1": {
                "torsions": {
                    "co-1-0": {"moltext": "moltext", "rel_energy": 10, "dihedral": "0"}
                },
                "delta_energy": 100,
            }
        }

        with self.subTest("key exists"):
            out = Template("{% load tags %}" "{{ data|get_delta_energy:1 }}").render(
                Context({"data": data})
            )
            self.assertEqual(out, "100.000")
        with self.subTest("key does not exist"):
            out = Template("{% load tags %}" "{{ data|get_delta_energy:2 }}").render(
                Context({"data": data})
            )
            self.assertEqual(out, "")

    def test_color_filter(self):
        """
        Test that the `color_filter` tag returns the expected color
        """
        with self.subTest("key exists"):
            data = {"HLM Prediction": "100"}
            out = Template(
                "{% load tags %}"
                "{% for key, value in data.items %}{{ value|color_filter:key }}{% endfor %}"
            ).render(Context({"data": data}))
            self.assertEqual(out, "red")
        with self.subTest("no assay_value but key exists"):
            data = {"assay_value": ""}
            out = Template(
                "{% load tags %}"
                "{% for key, value in data.items %}{{ value|color_filter:key }}{% endfor %}"
            ).render(Context({"data": data}))
            self.assertEqual(out, "None")

    def test_compared_color(self):
        """
        Test that the `compared_color` tag returns the expected color
        """
        with self.subTest("test red"):
            mmp = {"measured_data": {"assay_result": "10"}}
            compound = {"measured_data": {"assay_result": "8"}}
            out = Template(
                "{% load tags %}"
                "{{ mmp.measured_data.assay_result|compared_color:compound.measured_data.assay_result }}"
            ).render(Context({"mmp": mmp, "compound": compound}))
            self.assertEqual(out, "red")
        with self.subTest("test green"):
            mmp = {"measured_data": {"assay_result": "6"}}
            compound = {"measured_data": {"assay_result": "8"}}
            out = Template(
                "{% load tags %}"
                "{{ mmp.measured_data.assay_result|compared_color:compound.measured_data.assay_result }}"
            ).render(Context({"mmp": mmp, "compound": compound}))
            self.assertEqual(out, "green")
