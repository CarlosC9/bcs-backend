from biobarcoding.db_models.views_dashboards import View, DashboardView, Dashboard
from biobarcoding.rest import make_simple_rest_crud

bp_viewz, VizViewsAPI = make_simple_rest_crud(View, "viewz")
bp_dashboards, VizDashboardsAPI = make_simple_rest_crud(Dashboard, "dashboards")
bp_dashboard_items, VizDashboardItemsAPI = make_simple_rest_crud(DashboardView, "dashboard_items")
