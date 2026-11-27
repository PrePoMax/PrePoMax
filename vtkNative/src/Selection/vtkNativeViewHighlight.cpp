// PrePoMax - Copyright (C) 2016-2026 Matej Borovinsek
//
// The selection layer of vtkControl: the highlights of the mouse selection (the PickBy* hover paths and the
// highlight branch of PickByArea) and of the selected items, drawn by a renderer over the parts. A highlight is
// built from the data the Controller returns, as vtkMaxActor builds it (CreateNodalActor, CreatePolyFromData,
// CreateFromDataByVisualizationSurfaceExtraction), and drawn with the property vtkControl gives it.
#include "vtkNativeView.h"
#include "vtkNativeError.h"

#include <vtkCamera.h>
#include <vtkCellArray.h>
#include <vtkCellType.h>
#include <vtkDataSetMapper.h>
#include <vtkExtractEdges.h>
#include <vtkIdList.h>
#include <vtkNew.h>
#include <vtkPlaneCollection.h>
#include <vtkPointSet.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataSilhouette.h>
#include <vtkProperty.h>
#include <vtkTrivialProducer.h>
#include <vtkUnstructuredGridGeometryFilter.h>

#include <algorithm>
#include <vector>

namespace vtkn
{
    namespace
    {
        bool IsVolumeCell(int cellType)
        {
            return cellType == VTK_TETRA || cellType == VTK_HEXAHEDRON || cellType == VTK_WEDGE ||
                   cellType == VTK_PYRAMID || cellType == VTK_QUADRATIC_TETRA ||
                   cellType == VTK_QUADRATIC_HEXAHEDRON || cellType == VTK_QUADRATIC_WEDGE ||
                   cellType == VTK_QUADRATIC_PYRAMID;
        }

        vtkSmartPointer<vtkActor> CreateActor(vtkDataSet* data)
        {
            vtkNew<vtkDataSetMapper> mapper;
            mapper->SetInputData(data);
            mapper->ScalarVisibilityOff();
            vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
            actor->SetMapper(mapper);
            actor->PickableOff();
            return actor;
        }

        // VTK 9 lights a line or a point as if it faced the light and the camera, and draws an unlit one in the sum
        // of its ambient and diffuse colors, so the orange highlight of a line came out yellow; VTK 5.8 drew it in
        // about its own color, which an unlit actor with the color as its ambient color only gives
        bool HasOnlyLinesAndPoints(vtkDataSet* data)
        {
            for (vtkIdType i = 0; i < data->GetNumberOfCells(); i++)
            {
                int cellType = data->GetCellType(i);
                if (cellType != VTK_VERTEX && cellType != VTK_POLY_VERTEX && cellType != VTK_LINE &&
                    cellType != VTK_POLY_LINE && cellType != VTK_QUADRATIC_EDGE) return false;
            }
            return true;
        }

        // A vertex on every point a line of the data passes through
        vtkSmartPointer<vtkPolyData> GetLineJoints(vtkDataSet* data)
        {
            vtkSmartPointer<vtkPolyData> joints = vtkSmartPointer<vtkPolyData>::New();
            vtkPointSet* pointSet = vtkPointSet::SafeDownCast(data);
            if (pointSet == nullptr || pointSet->GetPoints() == nullptr) return joints;
            std::vector<bool> onLine(static_cast<size_t>(data->GetNumberOfPoints()), false);
            vtkNew<vtkIdList> pointIds;
            for (vtkIdType i = 0; i < data->GetNumberOfCells(); i++)
            {
                int cellType = data->GetCellType(i);
                if (cellType != VTK_LINE && cellType != VTK_POLY_LINE && cellType != VTK_QUADRATIC_EDGE) continue;
                data->GetCellPoints(i, pointIds);
                for (vtkIdType j = 0; j < pointIds->GetNumberOfIds(); j++) onLine[pointIds->GetId(j)] = true;
            }
            vtkNew<vtkCellArray> vertices;
            for (vtkIdType i = 0; i < static_cast<vtkIdType>(onLine.size()); i++)
                if (onLine[i]) vertices->InsertNextCell(1, &i);
            joints->SetPoints(pointSet->GetPoints());
            joints->SetVerts(vertices);
            return joints;
        }

        // Lines and round points in a flat color
        vtkSmartPointer<vtkActor> CreateUnlitActor(vtkDataSet* data, const double* color, double pointSize,
                                                   double lineWidth)
        {
            vtkSmartPointer<vtkActor> actor = CreateActor(data);
            vtkProperty* property = actor->GetProperty();
            property->SetColor(color[0], color[1], color[2]);
            property->SetLighting(false);
            property->SetAmbient(1);
            property->SetDiffuse(0);
            property->SetSpecular(0);
            property->SetPointSize(static_cast<float>(pointSize));
            property->SetLineWidth(static_cast<float>(lineWidth));
            property->BackfaceCullingOff();
            property->RenderPointsAsSpheresOn();
            DrawPointsAsFlatRoundDots(actor);
            return actor;
        }

        // How much wider than the lines and points the halo is, and how close to white its shade of their color
        const double HaloWidth = 2;
        const double HaloWhiteness = 0.65;

        void GetHaloColor(const double* color, double* haloColor)
        {
            for (int i = 0; i < 3; i++) haloColor[i] = color[i] + (1 - color[i]) * HaloWhiteness;
        }
    }

    // The selected lines drawn as the lines of the mouse highlight: with round joints and ends, and with the halo
    // under them. The added actors are cut by the clipping planes of the lines and follow their visibility. A symbol
    // is left as it is: its pipeline builds its lines again as the view zooms, and joints taken from the lines it had
    // would stay behind where they were
    void View::AddSelectionHalo(PartActor& actor)
    {
        const vtkn_part_style& style = actor.GetStyle();
        if (style.draw_on_geometry || style.keep_own_color) return;
        if (vtkTrivialProducer::SafeDownCast(actor.GetGeometry()->GetMapper()->GetInputAlgorithm()) == nullptr) return;
        vtkPolyData* polyData = actor.GetPolyData();
        if (polyData == nullptr || polyData->GetNumberOfLines() == 0 || polyData->GetNumberOfPolys() > 0 ||
            polyData->GetNumberOfStrips() > 0 || polyData->GetNumberOfVerts() > 0) return;
        vtkActor* geometry = actor.GetGeometry();
        vtkMapper* mapper = geometry->GetMapper();
        if (mapper->GetClippingPlanes() == nullptr)
        {
            vtkNew<vtkPlaneCollection> planes;
            mapper->SetClippingPlanes(planes);
        }
        double color[3];
        geometry->GetProperty()->GetColor(color);
        double lineWidth = geometry->GetProperty()->GetLineWidth();
        double haloColor[3];
        GetHaloColor(color, haloColor);
        double haloLineWidth = lineWidth + HaloWidth;
        vtkSmartPointer<vtkPolyData> lineJoints = GetLineJoints(polyData);
        AddSelectionHaloActor(CreateUnlitActor(polyData, haloColor, haloLineWidth, haloLineWidth), HaloRenderer,
                              geometry);
        if (lineJoints->GetNumberOfVerts() == 0) return;
        AddSelectionHaloActor(CreateUnlitActor(lineJoints, haloColor, haloLineWidth, haloLineWidth), HaloRenderer,
                              geometry);
        AddSelectionHaloActor(CreateUnlitActor(lineJoints, color, lineWidth, lineWidth), GetLayerRenderer(actor),
                              geometry);
    }

    void View::AddSelectionHaloActor(vtkActor* halo, vtkRenderer* renderer, vtkActor* source)
    {
        halo->PickableOff();
        halo->GetMapper()->SetClippingPlanes(source->GetMapper()->GetClippingPlanes());
        halo->SetPosition(source->GetPosition());
        halo->SetOrientation(source->GetOrientation());
        halo->SetScale(source->GetScale());
        halo->SetUserTransform(source->GetUserTransform());
        halo->SetVisibility(source->GetVisibility());
        renderer->AddActor(halo);
        SelectionHalos.push_back({ halo, renderer, source });
    }

    void View::ClearSelectionHalos()
    {
        for (const SelectionHaloEntry& entry : SelectionHalos) entry.Renderer->RemoveActor(entry.Actor);
        SelectionHalos.clear();
    }

    // Before the halo layer is drawn, the halos of the selected lines take the visibility and the position of their
    // lines
    void View::OnHaloRendererStart(vtkObject*, unsigned long, void* clientData, void*)
    {
        View* view = static_cast<View*>(clientData);
        for (const SelectionHaloEntry& entry : view->SelectionHalos)
        {
            if (entry.Actor->GetVisibility() != entry.Source->GetVisibility())
                entry.Actor->SetVisibility(entry.Source->GetVisibility());
            entry.Actor->SetPosition(entry.Source->GetPosition());
        }
    }

    // The selection renderer of vtkControl, the layer between the overlay and the widgets (InitializeGeometry)
    vtkRenderer* View::GetHighlightRenderer()
    {
        return SelectionRenderer;
    }

    // ApplySelectionFormattingToActor and the property the PickBy* paths set, as the caller computed them
    void View::AddHighlightActor(vtkActor* actor, const vtkn_highlight_style& style)
    {
        vtkProperty* property = actor->GetProperty();
        property->SetColor(style.color[0], style.color[1], style.color[2]);
        property->SetAmbient(style.ambient);
        property->SetDiffuse(style.diffuse);
        property->SetSpecular(style.specular);
        property->SetPointSize(static_cast<float>(style.point_size));
        property->SetLineWidth(static_cast<float>(style.line_width));
        property->SetOpacity(1);
        property->BackfaceCullingOff();
        vtkDataSet* data = actor->GetMapper()->GetInput();
        bool linesAndPoints = HasOnlyLinesAndPoints(data);
        if (linesAndPoints)
        {
            property->SetLighting(false);
            property->SetAmbient(1);
            property->SetDiffuse(0);
            property->SetSpecular(0);
        }
        vtkRenderer* renderer = style.on_geometry ? Renderer.GetPointer() : GetHighlightRenderer();
        if (!linesAndPoints || style.on_geometry)
        {
            renderer->AddActor(actor);
            Highlights.push_back({ actor, renderer, style.current != 0 });
            return;
        }
        // Round points, and a round dot as wide as the lines on every point of a line, which closes the corners
        // of thick lines and rounds their ends
        property->RenderPointsAsSpheresOn();
        DrawPointsAsFlatRoundDots(actor);
        renderer->AddActor(actor);
        Highlights.push_back({ actor, renderer, style.current != 0 });
        vtkSmartPointer<vtkPolyData> lineJoints = GetLineJoints(data);
        if (lineJoints->GetNumberOfVerts() > 0)
        {
            vtkSmartPointer<vtkActor> joints = CreateUnlitActor(lineJoints, style.color, style.line_width,
                                                                style.line_width);
            renderer->AddActor(joints);
            Highlights.push_back({ joints, renderer, style.current != 0 });
        }
        if (style.current != 0) AddHighlightHalo(data, lineJoints, style);
    }

    // The lines and points of the mouse highlight wider in a pale shade of its color, on the layer under it, which
    // sets them off from faces of any color
    void View::AddHighlightHalo(vtkDataSet* data, vtkPolyData* lineJoints, const vtkn_highlight_style& style)
    {
        double color[3];
        GetHaloColor(style.color, color);
        double lineWidth = style.line_width + HaloWidth;
        vtkSmartPointer<vtkActor> halo = CreateUnlitActor(data, color, style.point_size + HaloWidth, lineWidth);
        HaloRenderer->AddActor(halo);
        Highlights.push_back({ halo, HaloRenderer.GetPointer(), true });
        if (lineJoints->GetNumberOfVerts() == 0) return;
        vtkSmartPointer<vtkActor> joints = CreateUnlitActor(lineJoints, color, lineWidth, lineWidth);
        HaloRenderer->AddActor(joints);
        Highlights.push_back({ joints, HaloRenderer.GetPointer(), true });
    }

    void View::AddHighlight(const vtkn_mesh& mesh, const vtkn_highlight_style& style)
    {
        CheckMesh(mesh, "highlight");
        vtkNew<vtkPoints> points;
        points->SetNumberOfPoints(mesh.point_count);
        for (int i = 0; i < mesh.point_count; i++) points->SetPoint(i, mesh.coordinates + 3 * i);
        vtkNew<vtkUnstructuredGrid> grid;
        grid->SetPoints(points);
        bool hasVolumeCells = false;
        if (style.points || mesh.cell_count == 0)
        {
            // CreateNodalActor: a vertex per node
            grid->Allocate(mesh.point_count);
            for (vtkIdType i = 0; i < mesh.point_count; i++) grid->InsertNextCell(VTK_VERTEX, 1, &i);
        }
        else
        {
            grid->Allocate(mesh.cell_count);
            std::vector<vtkIdType> ids;
            for (int i = 0, cell = 0; i < mesh.cells_length; cell++)
            {
                int count = mesh.cells[i];
                ids.assign(mesh.cells + i + 1, mesh.cells + i + 1 + count);
                grid->InsertNextCell(mesh.cell_types[cell], count, ids.data());
                hasVolumeCells = hasVolumeCells || IsVolumeCell(mesh.cell_types[cell]);
                i += count + 1;
            }
        }
        // The visualization surface of volume cells: their outer faces
        vtkSmartPointer<vtkDataSet> surface = grid.GetPointer();
        if (hasVolumeCells)
        {
            vtkNew<vtkUnstructuredGridGeometryFilter> filter;
            filter->SetInputData(grid);
            filter->Update();
            surface = filter->GetOutput();
        }
        vtkSmartPointer<vtkActor> actor = CreateActor(surface);
        AddHighlightActor(actor, style);
        // The element edges keep the look of the part's element edges and are drawn on the same layer
        if (style.element_edges && !style.points)
        {
            vtkNew<vtkExtractEdges> edges;
            edges->SetInputData(surface);
            edges->Update();
            vtkSmartPointer<vtkActor> edgesActor = CreateActor(edges->GetOutput());
            vtkProperty* property = edgesActor->GetProperty();
            property->SetColor(0, 0, 0);
            property->SetLighting(false);
            property->SetLineWidth(1);
            vtkRenderer* renderer = style.on_geometry ? Renderer.GetPointer() : GetHighlightRenderer();
            renderer->AddActor(edgesActor);
            Highlights.push_back({ edgesActor, renderer, style.current != 0 });
        }
    }

    // GetHighlightActorFromActorName: a copy of the part's model edges, or its silhouette from the current camera
    void View::AddPartHighlight(const std::string& name, const vtkn_highlight_style& style)
    {
        auto found = Parts.find(name);
        if (found == Parts.end()) throw Error(VTKN_ERROR_NOT_FOUND, "No part is named '" + name + "'.");
        PartActor* part = found->second.get();
        vtkDataSet* modelEdges = part->GetModelEdges() != nullptr ? part->GetModelEdges()->GetMapper()->GetInput()
                                                                  : nullptr;
        vtkSmartPointer<vtkActor> actor;
        if (modelEdges != nullptr && modelEdges->GetNumberOfPoints() > 0) actor = CreateActor(modelEdges);
        else
        {
            vtkPolyData* surface = vtkPolyData::SafeDownCast(part->GetGeometry()->GetMapper()->GetInput());
            Require(surface != nullptr, "The part has no surface.");
            vtkNew<vtkPolyDataSilhouette> silhouette;
            silhouette->SetInputData(surface);
            silhouette->SetCamera(Renderer->GetActiveCamera());
            silhouette->SetDirectionToCameraVector();
            silhouette->SetBorderEdges(1);
            silhouette->SetEnableFeatureAngle(1);
            silhouette->SetFeatureAngle(30);
            silhouette->Update();
            vtkNew<vtkPolyData> copy;
            copy->DeepCopy(silhouette->GetOutput());
            actor = CreateActor(copy);
        }
        AddHighlightActor(actor, style);
    }

    void View::ClearHighlights(bool currentOnly)
    {
        auto removed = std::remove_if(Highlights.begin(), Highlights.end(), [currentOnly](const HighlightEntry& highlight) {
            if (currentOnly && !highlight.Current) return false;
            highlight.Renderer->RemoveActor(highlight.Actor);
            return true;
        });
        Highlights.erase(removed, Highlights.end());
    }
}
